% This function solves SDPs with unit-diagonal and without extra affine constraints using the primal approach:
% Min  <C, X>
% s.t. X >= 0
%      diag(X) = 1

function [X, obj, data] = ManiSDP_onlyunitdiag(C, options)

if nargin < 2 || isempty(options)
    options = struct();
end

if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 20; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-1; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.5; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 100; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 40; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
n = size(C,1);
fprintf('SDP size: n = %i, m = %i\n', n, n);

p = options.p0;
tau = options.tau;
Y = [];
U = [];
YC = [];
eG = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm; % tolerance on gradient norm

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    z = sum(C.*X);
    obj = full(sum(z));
    S = C - diag(z);
    [vS, dS] = eig(full(S), 'vector');
    dinf = max(0, -dS(1))/(1+dS(end));
    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    r = sum(e >= options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, dinf:%0.1e, r:%d, p:%d, time:%0.2fs\n', ...
             iter,    obj,       dinf,       r,    p,    toc(timespend));
    if dinf < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 20) == 0
        if iter > 50 && dinf > dinf0
            data.status = 2;
            fprintf('Slow progress!\n');
            break;
        else
            dinf0 = dinf;
        end
    end
    if r <= p - 1         
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    if p == 1
        nne = 1;
    else
        tau = max(0.4*tau, 1e-8);
        nne = min(sum(dS < -tau), options.delta);
    end
    if options.line_search == 1
        U = [zeros(p, n); vS(:,1:nne)'];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y; zeros(nne,n)];
    else
        Y = [Y; options.alpha*vS(:,1:nne)'];
        Y = Y./sqrt(sum(Y.^2));
    end
end
data.X = X;
data.S = S;
data.z = z;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
if data.status == 0 && dinf > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    function val = co(Y)
        val = sum((Y*C).*Y,'all');
    end

    function nY = line_search(Y, U)
        % The staircase direction satisfies Y'*U = 0, so the first
        % derivative vanishes. Use a second-order Armijo condition for
        % co(Retr_Y(alpha*U)), whose curvature is twice that of cost.
        nY = Y;
        % Rank truncation can perturb the unit-column constraints.
        ls_ynorm = sqrt(sum(Y.^2, 1));
        if any(~isfinite(ls_ynorm)) || any(ls_ynorm == 0)
            return;
        end
        Y = Y./ls_ynorm;
        nY = Y;
        ls_cost0 = co(Y);
        ls_u2 = sum(U.^2, 1);
        ls_columncost = sum((Y*C).*Y, 1);
        ls_curvature = 2*(sum((U*C).*U, 'all') ...
                         - sum(ls_columncost.*ls_u2));
        if ~isfinite(ls_cost0) || ~isfinite(ls_curvature) ...
                || ls_curvature >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(max(ls_u2)));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_trial = Y + ls_alpha*U;
            ls_trial = ls_trial./sqrt(sum(ls_trial.^2, 1));
            ls_cost = co(ls_trial);
            if isfinite(ls_cost) && ls_cost < ls_cost0 ...
                    && ls_cost - ls_cost0 <= ...
                       0.5*ls_armijo*ls_alpha^2*ls_curvature
                nY = ls_trial;
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep Y if no trial satisfies sufficient decrease.
    end
    
    function [f, store] = cost(Y, store)
        YC = Y*C;
        eG = sum(YC.*Y);
        f = 0.5*sum(eG);
    end

    function [G, store] = grad(Y, store)
        G = YC - Y.*eG;
    end

    function [H, store] = hess(Y, U, store)
        eH = U*C;
        H = eH - Y.*sum(Y.*eH) - U.*eG;
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) sum(d1.*d2, 'all');
        M.norm = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U);

        M.retr = @retraction;
        function y = retraction(x, d)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2));
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d);

        function x = random(n, m)
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end