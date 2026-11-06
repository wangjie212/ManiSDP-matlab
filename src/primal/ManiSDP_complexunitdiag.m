% Real SDP with paired unit diagonal, including data from convertCtoR:
% min <C,X>, A(X)=b, X>=0, X(j,j)+X(n+j,n+j)=1, j=1,...,n.
function [X, obj, data] = ManiSDP_complexunitdiag(At, b, c, K, options)
% Defaults selected on example_bqp_complex with d=30,40,50, tol=1e-8:
% p0=2, sigma0=sigma_min=1e-2, gamma=2, theta=1e-3, delta=8,
% TR_maxiter=4, TR_maxinner=25, kappa1=kappa2=1, line_search=1.
% The output X and dual certificate retain the original real SDP units.
% Exact complex blocks are handled internally in half the matrix dimension;
% general paired-diagonal data use the full real representation.
if nargin < 5 || isempty(options), options = struct(); end
if ~isstruct(K) || ~isfield(K,'s') || ~isscalar(K.s) ...
        || ~isreal(K.s) || ~isfinite(K.s) || K.s < 2 || mod(K.s,2) ~= 0
    error('K.s must be a positive even integer.');
end
n2 = K.s; n = n2/2;
b = b(:); c = c(:);
if size(At,1) ~= n2*n2 || size(At,2) ~= numel(b) || numel(c) ~= n2*n2 ...
        || ~isreal(At) || ~isreal(b) || ~isreal(c)
    error('Expected real SeDuMi data with size(At) = [K.s^2, numel(b)].');
end
if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 500; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-2; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 25; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1; end
if ~isfield(options,'line_search'); options.line_search = 1; end
if ~isfield(options,'verbosity'); options.verbosity = 1; end

timespend = tic;
At_original = At; c_original = c;
% Only symmetric parts contribute. Prepare them once, rather than
% symmetrizing a dense matrix at every Hessian application.
transpose_index = reshape(1:n2*n2,n2,n2).';
At = 0.5*(At + At(transpose_index(:),:));
C = reshape(c,n2,n2); C = (C+C')/2;
% For exact complex blocks, Z = Y1-i*Y2 is an isometry of the paired
% manifold. <[R,-I;I,R],Y'*Y> = real(<R+i*I,Z'*Z>). The real dual
% eigenvalues equal the complex dual eigenvalues, each repeated twice.
i11 = (1:n).' + n2*(0:n-1);
At11 = At(i11(:),:); At21 = At(i11(:)+n,:);
block_complex = isequal(At(i11(:)+n+n*n2,:),At11) ...
    && isequal(At(i11(:)+n*n2,:),-At21) ...
    && isequal(C(n+1:n2,n+1:n2),C(1:n,1:n)) ...
    && isequal(C(1:n,n+1:n2),-C(n+1:n2,1:n));
if block_complex
    At = At11 + 1i*At21;
    C = C(1:n,1:n) + 1i*C(n+1:n2,1:n);
    q = n;
else
    q = n2;
end
c = C(:); A = At';
% Release setup-only arrays before allocating trust-region workspaces.
C = []; At11 = []; At21 = []; i11 = []; transpose_index = []; %#ok<NASGU>
if options.verbosity
    fprintf('ManiSDP is starting...\n');
    fprintf('SDP size: n = %d, m = %d, internal matrix size = %d\n',n2,numel(b),q);
end
warning_state = warning('off','manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));
p = min(q,max(1,round(options.p0)));
sigma = min(options.sigma_max,max(realmin,options.sigma0));
tau = options.tau; y = zeros(numel(b),1); normb = 1+norm(b);
Y = []; fac_size = zeros(options.AL_maxiter,1);
diag_work = (1:q+1:q*q).';
grad_floor = min(options.tolgradnorm,0.1*options.tol);
previous_pinf = inf;
inner_iter = options.TR_maxiter; inner_calls = options.TR_maxinner;
problem.cost = @cost; problem.grad = @grad; problem.hess = @hess;
opts.verbosity = 0;
data.status = 0;
for iter = 1:options.AL_maxiter
    fac_size(iter) = p;
    problem.M = manifold(p);
    dual_c = c-At*y;
    opts.tolgradnorm = max(grad_floor,min(tau,0.1*previous_pinf));
    opts.maxiter = inner_iter;
    opts.maxinner = min(inner_calls,problem.M.dim());
    [Y,~,info] = trustregions(problem,Y,opts);
    gradnorm = info(end).gradnorm;
    gram = Y'*Y;
    obj = real(c'*gram(:));
    residual = real(A*gram(:))-b;
    pinf = norm(residual)/normb;
    pairinf = max(abs(point_inner(Y,Y)-1));
    y = y-sigma*residual;
    eS = reshape(c-At*y,q,q); eS = (eS+eS')/2;
    z_diag = real(sum(conj(gram).*eS,1));
    if block_complex, z_pair = z_diag;
    else, z_pair = z_diag(1:n)+z_diag(n+1:n2); end
    S_work = eS;
    if block_complex, z_work = z_pair;
    else, z_work = [z_pair,z_pair]; end
    S_work(diag_work) = S_work(diag_work)-z_work.';
    if ~isfinite(obj) || ~isfinite(pinf) || ~isfinite(gradnorm) ...
            || any(~isfinite(S_work(:)))
        dinf = inf; data.status = 2;
        break;
    end
    [vS,dS] = eig(S_work,'vector');
    [dS,order] = sort(real(dS)); vS = vS(:,order);
    dinf = max(0,-dS(1))/(1+max(0,dS(end)));
    by = b'*y+sum(z_pair);
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    eta = max([gap,pinf,dinf,pairinf]);
    [~,D,V] = svd(Y,'econ'); e = diag(D);
    rank_tol = min([options.theta,1e-3,sqrt(max(options.tol,min(1,eta)))]);
    r = max(1,sum(e >= rank_tol*e(1)));
    if options.verbosity
        fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
            iter,obj,gap,pinf,dinf,gradnorm,r,p,sigma,toc(timespend));
    end
    if eta < options.tol
        % Verify the actual real primal before accepting a smaller
        % internal representation's certificate.
        X_trial = real_primal(Y);
        trial_obj = c_original'*X_trial(:);
        trial_pinf = norm(At_original'*X_trial(:)-b)/normb;
        trial_gap = abs(trial_obj-by)/(1+abs(trial_obj)+abs(by));
        trial_diag = diag(X_trial);
        trial_pair = max(abs(trial_diag(1:n)+trial_diag(n+1:n2)-1));
        if max([trial_gap,trial_pinf,dinf,trial_pair]) < options.tol
            if options.verbosity, fprintf('Optimality is reached!\n'); end
            break;
        end
    end
    % Return the last solved point on exhaustion, without taking an
    % unverified extra rank step after the final certificate.
    if iter == options.AL_maxiter, break; end
    inner_inaccurate = gradnorm > max(opts.tolgradnorm,0.1*norm(residual));
    if inner_inaccurate && pinf > 0.9*previous_pinf && info(end).iter >= opts.maxiter
        inner_iter = min(max(4,options.TR_maxiter),2*inner_iter);
        inner_calls = min(max(30,options.TR_maxinner),2*inner_calls);
    elseif gradnorm <= opts.tolgradnorm || pinf < 0.7*previous_pinf
        inner_iter = max(options.TR_maxiter,ceil(inner_iter/2));
        inner_calls = max(options.TR_maxinner,ceil(inner_calls/2));
    end
    if r < p
        Y = normalize_point(V(:,1:r)'.*e(1:r));
        p = r;
    end
    tau = max(0.4*tau,grad_floor);
    nne = min([sum(dS < -tau),options.delta,q-p]);
    if nne > 0
        U = [zeros(p,q);vS(:,1:nne)'];
        Y = [Y;zeros(nne,q)]; %#ok<AGROW> Staircase changes the factor size.
        % dual_c still represents the just-solved augmented objective.
        Y = line_search(Y,U);
        p = p+nne;
    end
    if pinf < options.kappa1*gradnorm
        sigma = max(sigma/options.gamma,options.sigma_min);
    elseif pinf > options.kappa2*gradnorm
        sigma = min(sigma*options.gamma,options.sigma_max);
    end
    previous_pinf = pinf;
end
X = real_primal(Y);
if block_complex
    S = [real(S_work),-imag(S_work);imag(S_work),real(S_work)];
else
    S = real(S_work);
end
obj = c_original'*X(:);
pinf = norm(At_original'*X(:)-b)/normb;
by = b'*y+sum(z_pair);
gap = abs(obj-by)/(1+abs(obj)+abs(by));
diagonal = diag(X);
pairinf = max(abs(diagonal(1:n)+diagonal(n+1:n2)-1));
eta = max([gap,pinf,dinf,pairinf]);
if data.status == 0 && eta >= options.tol, data.status = 1; end
data.X = X; data.y = y; data.S = S; data.z = [z_pair,z_pair];
data.gap = gap; data.pinf = pinf; data.dinf = dinf;
data.pairinf = pairinf; data.gradnorm = gradnorm;
data.fac_size = fac_size(1:iter); data.block_complex = block_complex;
data.time = toc(timespend);
if options.verbosity
    if data.status == 1, fprintf('Iteration maximum is reached!\n'); end
    fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n',obj,data.time);
end

    function Xr = real_primal(point)
        if block_complex, point = [real(point),-imag(point)]; end
        Xr = point'*point;
    end

    function lam = point_inner(point,direction)
        if block_complex
            lam = real(sum(conj(point).*direction,1));
        else
            paired = sum(point(:,1:n).*direction(:,1:n),1) ...
                + sum(point(:,n+1:n2).*direction(:,n+1:n2),1);
            lam = [paired,paired];
        end
    end

    function point = normalize_point(point)
        point = point./sqrt(point_inner(point,point));
    end

    function tangent = project(point,direction)
        tangent = direction-point.*point_inner(point,direction);
    end

    function M = manifold(rank)
        M.dim = @() (2*rank-1)*n;
        M.inner = @(x,u,v) real(u(:)'*v(:));
        M.norm = @(x,u) norm(u,'fro');
        M.typicaldist = @() pi*sqrt(n);
        M.proj = @project; M.tangent = @project;
        M.retr = @(point,direction) normalize_point(point+direction);
        M.rand = @() random_point(rank);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(rank,q);
        M.transp = @(x1,x2,u) project(x2,u);
    end

    function point = random_point(rank)
        if block_complex, point = randn(rank,n)+1i*randn(rank,n);
        else, point = [randn(rank,n),randn(rank,n)]; end
        point = normalize_point(point);
    end

    function [f,store] = cost(point,store)
        g = point'*point;
        store.residual = real(A*g(:))-b;
        % Multiplier-dependent constants need not enter cost differences.
        f = real(dual_c'*g(:))+0.5*sigma*(store.residual'*store.residual);
    end

    function [G,store] = grad(point,store)
        if ~isfield(store,'residual'), [~,store] = cost(point,store); end
        store.slack = reshape(dual_c+sigma*(At*store.residual),q,q);
        eG = 2*point*store.slack;
        store.lam = point_inner(point,eG);
        G = eG-point.*store.lam;
    end

    function [H,store] = hess(point,direction,store)
        if ~isfield(store,'slack') || ~isfield(store,'lam')
            [~,store] = grad(point,store);
        end
        yu = point'*direction;
        ayu = reshape(At*real(A*yu(:)),q,q);
        eH = 2*direction*store.slack+4*sigma*(point*ayu);
        H = project(point,eH)-direction.*store.lam;
    end

    function next = line_search(point,direction)
        point = normalize_point(point); next = point;
        u2 = point_inner(direction,direction);
        if any(~isfinite(u2)) || max(u2) == 0, return; end
        base = point'*point; utu = direction'*direction;
        base = correct_diagonal(base,true);
        dx = utu-0.5*(u2.'.*base+base.*u2);
        dx = correct_diagonal(dx,false);
        residual0 = real(A*base(:))-b;
        linear = dual_c+sigma*(At*residual0);
        if block_complex
            linear(diag_work) = 0;
        else
            first = diag_work(1:n); second = diag_work(n+1:n2);
            diff = 0.5*(linear(first)-linear(second));
            linear(first) = diff; linear(second) = -diff;
        end
        slope = real(linear'*dx(:));
        if ~isfinite(slope) || slope >= 0, return; end
        alpha = min(1,1/sqrt(max(u2)));
        if options.line_search == 0, alpha = min(alpha,options.alpha); end
        for ls_iter = 1:50
            theta = alpha^2; scale = sqrt(1+theta*u2);
            shift = theta*u2./(1+scale);
            dx = (theta*utu-(shift.'+shift+shift.'*shift).*base)./(scale.'*scale);
            dx = correct_diagonal(dx,false);
            adx = real(A*dx(:));
            change = real(linear'*dx(:))+0.5*sigma*(adx'*adx);
            if isfinite(change) && change < 0 && change <= 0.1*theta*slope
                next = normalize_point(point+alpha*direction);
                return;
            end
            alpha = 0.5*alpha;
        end
    end

    function g = correct_diagonal(g,basepoint)
        if block_complex
            g(diag_work) = double(basepoint);
        else
            first = diag_work(1:n); second = diag_work(n+1:n2);
            g(second) = double(basepoint)-g(first);
        end
    end
end

