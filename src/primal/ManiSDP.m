% Primal real SDP: min <C,X>, At'*X(:)=b, X positive semidefinite.
function [X,obj,data]=ManiSDP(At,b,c,K,options)
% Best tested defaults on example_qsphere, d=30,40,50 (tol=1e-8):
% p0=1, sigma0=sigma_min=1e-2, gamma=2, TR_maxiter=4,
% TR_maxinner=20, kappa1=0.02, kappa2=0.1, theta=1e-2, delta=8.
% Validated on problem seeds 1:10 per size, with solver seed 0.
% Both staircase modes enforce decrease of the current augmented cost;
% line_search=0 uses alpha as an initial bound. Large theta is capped.
% Inner budgets are starting values: stalled small budgets recover up to
% 4 outer trust-region steps and 20 inner CG steps before being restored.
% Derivative caches use Manopt point keys and reset for each subproblem.
% data.status: 0 = converged, 1 = iteration limit, 2 = nonfinite iterate.
if nargin<5||isempty(options),options=struct();end
if ~isstruct(K)||~isfield(K,'s')||~isscalar(K.s)||~isreal(K.s) ...
        ||~isfinite(K.s)||K.s<1||K.s~=fix(K.s)
    error('K.s must be a positive scalar integer.');
end
n=K.s;b=b(:);c=c(:);m=numel(b);
if size(At,1)~=n*n||size(At,2)~=m||numel(c)~=n*n ...
        ||~isreal(At)||~isreal(b)||~isreal(c)
    error('Expected real SDP data with size(At)=[K.s^2,numel(b)].');
end
if ~isfield(options,'p0'),options.p0=1;end
if ~isfield(options,'AL_maxiter'),options.AL_maxiter=1000;end
if ~isfield(options,'gamma'),options.gamma=2;end
if ~isfield(options,'tau'),options.tau=1e-2;end
if ~isfield(options,'sigma0'),options.sigma0=1e-2;end
if ~isfield(options,'sigma_min'),options.sigma_min=1e-2;end
if ~isfield(options,'sigma_max'),options.sigma_max=1e7;end
if ~isfield(options,'tol'),options.tol=1e-8;end
if ~isfield(options,'theta'),options.theta=1e-2;end
if ~isfield(options,'delta'),options.delta=8;end
if ~isfield(options,'alpha'),options.alpha=0.1;end
if ~isfield(options,'tolgradnorm'),options.tolgradnorm=1e-8;end
if ~isfield(options,'TR_maxinner'),options.TR_maxinner=20;end
if ~isfield(options,'TR_maxiter'),options.TR_maxiter=4;end
if ~isfield(options,'kappa1'),options.kappa1=0.02;end
if ~isfield(options,'kappa2'),options.kappa2=0.1;end
if ~isfield(options,'line_search'),options.line_search=1;end
if ~isfield(options,'verbosity'),options.verbosity=1;end
if options.AL_maxiter<1||options.AL_maxiter~=fix(options.AL_maxiter) ...
        ||options.sigma_min<=0||options.sigma_max<options.sigma_min ...
        ||options.sigma0<=0||options.gamma<=1||options.tol<=0 ...
        ||options.tolgradnorm<=0||options.TR_maxiter<1||options.TR_maxinner<1
    error('Budgets, tolerances and penalties must be positive, with gamma > 1.');
end
timespend=tic;
transpose_index=reshape(1:n*n,n,n).';
At=0.5*(At+At(transpose_index(:),:));A=At';
C=reshape(c,n,n);C=(C+C')/2;c=C(:);
normb=1+norm(b);p=min(n,max(1,round(options.p0)));
sigma=min(options.sigma_max,max(options.sigma_min,options.sigma0));
y=zeros(m,1);tau=options.tau;Y=[];U=[];
if isfield(options,'Y0')&&~isempty(options.Y0)
    Y=options.Y0;
    if size(Y,1)~=n||~isreal(Y)||any(~isfinite(Y(:)))
        error('Y0 must be a finite real factor with K.s rows.');
    end
    p=size(Y,2);
    if p>n,[v0,d0,~]=svd(Y,'econ');Y=v0*d0;p=n;end
end
fac_size=zeros(options.AL_maxiter,1);seta=fac_size;
grad_floor=min(options.tolgradnorm,options.tol);
scaled_At=[];scaled_penalty=NaN;
cached_slack=[];cached_key='';
previous_pinf=inf;inner_iter=options.TR_maxiter;inner_calls=options.TR_maxinner;
problem.cost=@cost;problem.grad=@grad;problem.hess=@hess;
opts.verbosity=0;data.status=0;
if options.verbosity,fprintf('ManiSDP is starting...\nSDP size: n = %d, m = %d\n',n,m);end
for iter=1:options.AL_maxiter
    fac_size(iter)=p;problem.M=euclideanfactory(n,p);
    cached_slack=[];cached_key='';S=[]; %#ok<NASGU> Release the prior dual matrix.
    if sigma~=scaled_penalty,scaled_At=sigma*At;scaled_penalty=sigma;end
    if ~isempty(U),Y=line_search(Y,U);end
    opts.tolgradnorm=max(grad_floor,tau);
    opts.maxiter=inner_iter;opts.maxinner=min(inner_calls,n*p);
    [Y,~,info]=trustregions(problem,Y,opts);
    cached_slack=[];
    gradnorm=info(end).gradnorm;
    X=Y*Y';x=X(:);residual=A*x-b;
    pinf=norm(residual)/normb;
    y=y-sigma*residual;
    obj=c'*x;by=b'*y;gap=abs(obj-by)/(1+abs(obj)+abs(by));
    S=reshape(c-At*y,n,n);S=(S+S')/2;
    if ~isfinite(obj)||~isfinite(pinf)||~isfinite(gradnorm)||any(~isfinite(S(:)))
        dinf=inf;seta(iter)=inf;data.status=2;break;
    end
    [vectors,values]=eig(S,'vector');
    [values,order]=sort(real(values));vectors=vectors(:,order);
    dinf=max(0,-values(1))/(1+max(0,values(end)));
    eta=max([pinf,gap,dinf]);seta(iter)=eta;
    [V,D,~]=svd(Y,'econ');e=diag(D);
    r=max(1,sum(e>=min(options.theta,1e-2)*e(1)));
    if options.verbosity
        fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3g, time:%0.2fs\n', ...
            iter,obj,gap,pinf,dinf,gradnorm,r,p,sigma,toc(timespend));
    end
    if eta<options.tol||iter==options.AL_maxiter,break;end
    % Recover useful budgets when an undersolved subproblem stalls.
    if (inner_iter<4||inner_calls<20)&&gradnorm>opts.tolgradnorm ...
            &&pinf>0.9*previous_pinf&&info(end).iter>=opts.maxiter
        inner_iter=min(max(4,options.TR_maxiter),max(4,2*inner_iter));
        inner_calls=min(max(20,options.TR_maxinner),max(20,2*inner_calls));
    elseif pinf<0.7*previous_pinf||gradnorm<=opts.tolgradnorm
        inner_iter=max(options.TR_maxiter,ceil(inner_iter/2));
        inner_calls=max(options.TR_maxinner,ceil(inner_calls/2));
    end
    if r<p,Y=V(:,1:r)*diag(e(1:r));p=r;end
    tau=max(0.4*tau,grad_floor);
    if norm(Y,'fro')==0&&values(1)<-tau,Y=zeros(n,0);p=0;end
    nne=min([sum(values<-tau),options.delta,n-p]);
    U=[zeros(n,p),vectors(:,1:nne)];
    Y=[Y,zeros(n,nne)]; %#ok<AGROW>
    p=p+nne;
    if pinf<options.kappa1*gradnorm
        sigma=max(sigma/options.gamma,options.sigma_min);
    elseif pinf>options.kappa2*gradnorm
        sigma=min(sigma*options.gamma,options.sigma_max);
    end
    previous_pinf=pinf;
end
eta=max([pinf,gap,dinf]);
if data.status==0&&eta>=options.tol,data.status=1;end
data.X=X;data.y=y;data.S=S;data.gap=gap;data.pinf=pinf;data.dinf=dinf;
data.gradnorm=gradnorm;data.time=toc(timespend);
data.fac_size=fac_size(1:iter);data.seta=seta(1:iter);
if options.verbosity
    fprintf('ManiSDP: status = %d, optimum = %0.8f, time = %0.2fs\n',data.status,obj,data.time);
end

    function [f,store]=cost(point,store)
        xx=point*point';
        store.residual=A*xx(:)-b-y/sigma;
        f=c'*xx(:)+sigma/2*(store.residual'*store.residual);
    end

    function G=grad(point,storedb,key)
        if ~strcmp(key,cached_key),cache_slack(point,storedb,key);end
        G=2*cached_slack*point;
    end

    function H=hess(point,direction,storedb,key)
        if ~strcmp(key,cached_key),cache_slack(point,storedb,key);end
        yu=direction*point';
        curvature=reshape(A'*(At'*yu(:)),n,n);
        H=2*cached_slack*direction+4*sigma*(curvature*point);
    end

    function cache_slack(point,storedb,key)
        store=storedb.getWithShared(key);
        if ~isfield(store,'residual')
            [~,store]=cost(point,store);
            storedb.setWithShared(store,key);
        end
        % A point identifier prevents stale derivatives after other trial
        % points are evaluated, without caching several dense matrices.
        cached_slack=reshape(c+scaled_At*store.residual,n,n);
        cached_key=key;
    end

    function next=line_search(point,direction)
        next=point;u2=sum(direction.^2,'all');
        if ~isfinite(u2)||u2==0||any(~isfinite(point(:))),return;end
        base=point*point';step_matrix=direction*direction';
        res=A*base(:)-b-y/sigma;Ad=A*step_matrix(:);
        slope=c'*step_matrix(:)+sigma*(res'*Ad);
        quadratic=0.5*sigma*(Ad'*Ad);
        if ~isfinite(slope)||slope>=0||~isfinite(quadratic),return;end
        alpha=min(1,1/sqrt(u2));
        if options.line_search==0,alpha=min(alpha,options.alpha);end
        for ls=1:50
            theta=alpha^2;change=theta*(slope+theta*quadratic);
            if isfinite(change)&&change<0&&change<=0.1*theta*slope
                next=point+alpha*direction;return;
            end
            alpha=alpha/2;
        end
    end
end
