% Dual unit-diagonal SDP:
% min b'*y, S = C + Apsd'*y >= 0, diag(S)=1, B'*y=cf,
% where the input matrix is [B,Apsd] and c=[cf;C(:)].
% The returned certificate is X=eX-Diag(z)>=0, Apsd*eX(:)+B*w=b.
% Its objective is -<C,eX>+cf'*w+sum(z). For SOS BQP data C=0.
function [X, obj, data] = ManiDSDP_unitdiag(A, b, c, K, options)
% Defaults selected on example_bqp_dual, d=40,50,60, at tol=1e-8:
% Validated on 10 instances per size (problem seeds 1:10, solver seed 0).
% p0=ceil(log(m)), sigma0=sigma_min=3e-4, TR_maxiter=4,
% TR_maxinner=20, kappa1=10, kappa2=100, delta=8, line_search=1.
% Inner budgets adapt on stagnation. Both staircase modes backtrack;
% line_search=0 uses alpha as the initial trial-step bound.
if nargin < 5 || isempty(options), options=struct(); end
if ~isstruct(K) || ~isfield(K,'s') || ~isscalar(K.s) ...
        || ~isreal(K.s) || ~isfinite(K.s) || K.s<1 || K.s~=fix(K.s)
    error('K.s must be a positive scalar integer.');
end
if ~isfield(K,'f') || isempty(K.f), K.f=0; end
if ~isscalar(K.f) || ~isreal(K.f) || ~isfinite(K.f) || K.f<0 || K.f~=fix(K.f)
    error('K.f must be a nonnegative scalar integer.');
end
n=K.s; nf=K.f; b=b(:); c=c(:); m=numel(b);
if size(A,1)~=m || size(A,2)~=nf+n*n || numel(c)~=nf+n*n ...
        || ~isreal(A) || ~isreal(b) || ~isreal(c)
    error('Expected real SeDuMi data with size(A)=[numel(b),K.f+K.s^2].');
end
if ~isfield(options,'p0'); options.p0=max(2,ceil(log(max(m,2)))); end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter=300; end
if ~isfield(options,'gamma'); options.gamma=2; end
if ~isfield(options,'tau'); options.tau=1e-2; end
if ~isfield(options,'sigma0'); options.sigma0=3e-4; end
if ~isfield(options,'sigma_min'); options.sigma_min=3e-4; end
if ~isfield(options,'sigma_max'); options.sigma_max=1e7; end
if ~isfield(options,'tol'); options.tol=1e-8; end
if ~isfield(options,'theta'); options.theta=1e-3; end
if ~isfield(options,'delta'); options.delta=8; end
if ~isfield(options,'alpha'); options.alpha=0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm=1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner=20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter=4; end
if ~isfield(options,'kappa1'); options.kappa1=10; end
if ~isfield(options,'kappa2'); options.kappa2=100; end
if ~isfield(options,'line_search'); options.line_search=1; end
if ~isfield(options,'verbosity'); options.verbosity=1; end

timespend=tic;
normc=1+norm(c); normb=1+norm(b);
B=A(:,1:nf); A=A(:,nf+1:end); cf=reshape(c(1:nf),nf,1); c=c(nf+1:end);
% Only symmetric parts act on PSD matrices. Canonicalize once so the
% projection and the Hessian use the same self-adjoint operator.
transpose_index=reshape(1:n*n,n,n).';
A=0.5*(A+A(:,transpose_index(:)));
C=reshape(c,n,n); C=(C+C')/2; c=C(:);
diag_index=(1:n+1:n*n).';
% Validate the row Gram instead of trusting a possibly stale dAAt hint.
gram=A*A'; gram=0.5*(gram+gram');
gdiag=full(diag(gram)); iA=[]; chol_factor=[];
diagonal_gram=nnz(gram)==m;
if m==0
    diagonal_gram=true;
elseif diagonal_gram
    if any(~isfinite(gdiag)) || any(gdiag<=0)
        error('The PSD constraint block must have full row rank.');
    end
    iA=(spdiags(1./gdiag,0,m,m)*A)';
else
    [chol_factor,flag]=chol(gram);
    if flag~=0, error('The PSD constraint block must have full row rank.'); end
end
hint_corrected=false;
if isfield(options,'dAAt')
    hint=options.dAAt(:);
    hint_corrected=~diagonal_gram || numel(hint)~=m;
    if ~hint_corrected
        hint_corrected=any(~isfinite(hint)) || ...
            any(abs(hint-gdiag)>100*eps*max(1,abs(gdiag)));
    end
end
bA=lift(b); iAB=lift(B);
free_diag=iAB(diag_index,:);
diagonal_free=nnz(iAB)==nnz(free_diag);
gram=[]; C=[]; transpose_index=[]; %#ok<NASGU>
if options.verbosity
    fprintf('ManiDSDP is starting...\n');
    fprintf('SDP size: n = %d, m = %d\n',n,m);
end
warning_state=warning('off','manopt:trs_tCG_cached:memory');
warning_cleanup=onCleanup(@() warning(warning_state));
p=min(n,max(1,round(options.p0)));
sigma=min(options.sigma_max,max(options.sigma_min,options.sigma0));
tau=options.tau; x=zeros(n*n,1); w=zeros(nf,1); Y=[];
fac_size=zeros(options.ADMM_maxiter,1); seta=fac_size;
grad_floor=min(options.tolgradnorm,0.1*options.tol);
previous_pinf=inf;
inner_iter=options.TR_maxiter; inner_calls=options.TR_maxinner;
problem.cost=@cost; problem.grad=@grad; problem.hess=@hess;
opts.verbosity=0; data.status=0;
for iter=1:options.ADMM_maxiter
    fac_size(iter)=p; problem.M=manifold(p);
    % Freeze the augmented objective until its rank step has been tested.
    model_x=x; model_w=w; penalty=sigma;
    model_base=bA+x-project(x);
    if diagonal_free
        model_base(diag_index)=model_base(diag_index)-free_diag*w;
    else
        model_base=model_base-iAB*w;
    end
    opts.tolgradnorm=max(grad_floor,min(tau,0.1*previous_pinf));
    opts.maxiter=inner_iter;
    opts.maxinner=min(inner_calls,max(1,problem.M.dim()));
    [Y,~,info]=trustregions(problem,Y,opts);
    gradnorm=info(end).gradnorm;
    S=Y'*Y; sc=S(:)-c;
    y=dual_map(sc);
    As=A'*y-sc; Af=B'*y-cf;
    pinf=(norm(As)+norm(Af))/normc;
    % This is the exact augmented-objective gradient in S. Its update
    % preserves A*eX+B*w=b, including general (non-diagonal) free blocks.
    eX=gradient_matrix(As,Af);
    x=eX(:)-bA; w=model_w-penalty*Af;
    z=sum(S.*eX,1).';
    X=eX; X(diag_index)=X(diag_index)-z;
    obj=-c'*eX(:)+cf'*w+sum(z); by=b'*y;
    primalinf=norm(A*eX(:)+B*w-b)/normb;
    diaginf=max(abs(diag(S)-1));
    if ~isfinite(obj) || ~isfinite(pinf) || ~isfinite(gradnorm) ...
            || any(~isfinite(X(:)))
        dinf=inf; data.status=2;
        break;
    end
    [vX,dX]=eig(X,'vector');
    [dX,order]=sort(real(dX)); vX=vX(:,order);
    dinf=max(0,-dX(1))/(1+max(0,dX(end)));
    gap=abs(obj-by)/(1+abs(obj)+abs(by));
    eta=max([gap,pinf,dinf,primalinf,diaginf]); seta(iter)=eta;
    [~,D,V]=svd(Y,'econ'); e=diag(D);
    rank_tol=min([options.theta,1e-3,sqrt(max(options.tol,min(1,eta)))]);
    r=max(1,sum(e>=rank_tol*e(1)));
    if options.verbosity
        fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
            iter,obj,gap,pinf,dinf,gradnorm,r,p,sigma,toc(timespend));
    end
    if eta<options.tol
        if options.verbosity, fprintf('Optimality is reached!\n'); end
        break;
    end
    if iter==options.ADMM_maxiter, break; end
    inner_inaccurate=gradnorm>opts.tolgradnorm;
    low_budget=inner_iter<4 || inner_calls<20;
    if inner_inaccurate && pinf>0.9*previous_pinf && info(end).iter>=opts.maxiter ...
            && (low_budget || gradnorm>0.1*(norm(As)+norm(Af)))
        % A tiny requested budget must not trap the penalty update at an
        % inaccurate stationary point. Recover the normal budget at once.
        inner_iter=min(max(4,options.TR_maxiter),max(4,2*inner_iter));
        inner_calls=min(max(30,options.TR_maxinner),max(20,2*inner_calls));
    elseif gradnorm<=opts.tolgradnorm || pinf<0.7*previous_pinf
        inner_iter=max(options.TR_maxiter,ceil(inner_iter/2));
        inner_calls=max(options.TR_maxinner,ceil(inner_calls/2));
    end
    if r<p
        Y=normalize_point(V(:,1:r)'.*e(1:r)); p=r;
    end
    tau=max(0.4*tau,grad_floor);
    nne=min([sum(dX < -tau),options.delta,n-p]);
    if nne>0
        U=[zeros(p,n);vX(:,1:nne)'];
        Y=[Y;zeros(nne,n)]; %#ok<AGROW> Staircase changes the factor size.
        Y=line_search(Y,U);
        p=p+nne;
    end
    if pinf<options.kappa1*gradnorm
        sigma=max(sigma/options.gamma,options.sigma_min);
    elseif pinf>options.kappa2*gradnorm
        sigma=min(sigma*options.gamma,options.sigma_max);
    end
    previous_pinf=pinf;
end
gap=abs(obj-b'*y)/(1+abs(obj)+abs(b'*y));
eta=max([gap,pinf,dinf,primalinf,diaginf]);
seta(iter)=eta;
if data.status==0 && eta>=options.tol, data.status=1; end
data.X=X; data.eX=eX; data.z=z; data.y=y; data.S=S; data.w=w;
data.gap=gap; data.pinf=pinf; data.dinf=dinf;
data.primalinf=primalinf; data.diaginf=diaginf; data.gradnorm=gradnorm;
data.fac_size=fac_size(1:iter); data.seta=seta(1:iter);
data.diagonal_gram=diagonal_gram; data.dAAt_corrected=hint_corrected;
data.time=toc(timespend);
if options.verbosity
    if data.status==1, fprintf('Iteration maximum is reached!\n'); end
    fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n',obj,data.time);
end

    function out=dual_map(s)
        if m==0
            out=zeros(0,size(s,2));
        elseif diagonal_gram
            out=iA'*s;
        else
            out=chol_factor\(chol_factor'\(A*s));
        end
    end

    function out=lift(rhs)
        if m==0
            out=zeros(n*n,size(rhs,2));
        elseif diagonal_gram
            out=iA*rhs;
        else
            out=A'*(chol_factor\(chol_factor'\rhs));
        end
    end

    function out=project(s)
        out=A'*dual_map(s);
    end

    function E=gradient_matrix(as,af)
        vector=model_base-penalty*as;
        if diagonal_free
            vector(diag_index)=vector(diag_index)+penalty*(free_diag*af);
        else
            vector=vector+penalty*(iAB*af);
        end
        E=reshape(vector,n,n); E=(E+E')/2;
    end

    function point=normalize_point(point)
        point=point./sqrt(sum(point.^2,1));
    end

    function tangent=tangent_projection(point,direction)
        tangent=direction-point.*sum(point.*direction,1);
    end

    function M=manifold(rank)
        M.dim=@() (rank-1)*n;
        M.inner=@(x,u,v) u(:)'*v(:);
        M.norm=@(x,u) norm(u,'fro');
        M.typicaldist=@() pi*sqrt(n);
        M.proj=@tangent_projection; M.tangent=@tangent_projection;
        M.retr=@(point,direction) normalize_point(point+direction);
        M.rand=@() normalize_point(randn(rank,n));
        M.lincomb=@matrixlincomb; M.zerovec=@(x) zeros(rank,n);
        M.transp=@(x1,x2,u) tangent_projection(x2,u);
    end

    function [f,store]=cost(point,store)
        gram_s=point'*point; sc0=gram_s(:)-c;
        y0=dual_map(sc0);
        store.As=A'*y0-sc0;
        store.Af=B'*y0-cf;
        % Dropping multiplier-only constants stabilizes cost differences.
        f=b'*y0+0.5*penalty*(store.As'*store.As+store.Af'*store.Af) ...
            -model_x'*store.As-model_w'*store.Af;
    end

    function [G,store]=grad(point,store)
        if ~isfield(store,'E')
            if ~isfield(store,'As'), [~,store]=cost(point,store); end
            store.E=gradient_matrix(store.As,store.Af);
            % The Hessian needs E, not another dense copy of the residual.
            store=rmfield(store,{'As','Af'});
        end
        eG=2*point*store.E;
        store.lam=sum(point.*eG,1);
        G=eG-point.*store.lam;
    end

    function [H,store]=hess(point,direction,store)
        if ~isfield(store,'E') || ~isfield(store,'lam')
            [~,store]=grad(point,store);
        end
        yu=point'*direction;
        pyu=reshape(project(yu(:)),n,n);
        eH=2*direction*store.E-4*penalty*(point*pyu) ...
            +2*penalty*((point*direction')*point+(point*point')*direction);
        if diagonal_free
            diagonal_yu=sum(point.*direction,1).';
            curvature=free_diag*(free_diag'*diagonal_yu);
            eH=eH+4*penalty*(point.*curvature.');
        else
            free_yu=reshape(iAB*(iAB'*yu(:)),n,n);
            eH=eH+4*penalty*(point*free_yu);
        end
        H=tangent_projection(point,eH)-direction.*store.lam;
    end

    function next=line_search(point,direction)
        point=normalize_point(point); next=point;
        u2=sum(direction.^2,1);
        if any(~isfinite(u2)) || max(u2)==0, return; end
        base=point'*point; base(diag_index)=1;
        utu=direction'*direction;
        ds=utu-0.5*(u2.'.*base+base.*u2); ds(diag_index)=0;
        sc0=base(:)-c; y0=dual_map(sc0);
        linear=gradient_matrix(A'*y0-sc0,B'*y0-cf);
        linear(diag_index)=0;
        slope=linear(:)'*ds(:);
        if ~isfinite(slope) || slope>=0, return; end
        alpha=min(1,1/sqrt(max(u2)));
        if options.line_search==0, alpha=min(alpha,options.alpha); end
        for ls_iter=1:50
            theta=alpha^2; scale=sqrt(1+theta*u2);
            shift=theta*u2./(1+scale);
            ds=(theta*utu-(shift.'+shift+shift.'*shift).*base)./(scale.'*scale);
            ds(diag_index)=0;
            das=project(ds(:))-ds(:);
            if diagonal_free, daf=zeros(nf,1);else,daf=iAB'*ds(:);end
            change=linear(:)'*ds(:)+0.5*penalty*(das'*das+daf'*daf);
            if isfinite(change) && change<0 && change<=0.1*theta*slope
                next=normalize_point(point+alpha*direction);
                return;
            end
            alpha=0.5*alpha;
        end
    end
end

