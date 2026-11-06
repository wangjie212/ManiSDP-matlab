% This function solves SDPs with unit-diagonal using the primal approach:
% Min  <C, X>
% s.t. A(X) = b
%      X >= 0
%      diag(X) = 1

function [X, obj, data] = ManiSDP_unitdiag(At, b, c, K, options)
% Defaults were selected on example_bqp with d=30,40 and tol=1e-8:
% sigma0=1e-3, gamma=2, TR_maxiter=4, TR_maxinner=20,
% kappa1=0.5, kappa2=2, theta=1e-3, delta=8, line_search=1.
% TR budgets adapt on stagnation. Both staircase modes use sufficient
% decrease; line_search=0 uses alpha as its initial trial step.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s) ...
        || ~isfinite(K.s) || K.s < 1 || K.s ~= fix(K.s)
    error('K.s must be a positive scalar integer.');
end

if ~isfield(options,'p0'); options.p0 = max(2, ceil(log(max(length(b), 2)))); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-3; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-3; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 0.5; end
if ~isfield(options,'kappa2'); options.kappa2 = 2; end
if ~isfield(options,'line_search'); options.line_search = 1; end
if ~isfield(options,'verbosity'); options.verbosity = 1; end

n = K.s;
if options.verbosity
    fprintf('ManiSDP is starting...\n');
    fprintf('SDP size: n = %i, m = %i\n', n, numel(b));
end
warning_state = warning('off', 'manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));

A = At';
p = min(n, max(2, round(options.p0)));
sigma = min(options.sigma_max, max(options.sigma_min, options.sigma0));
gamma = options.gamma;
tau = options.tau;
y = zeros(length(b), 1);
normb = 1 + norm(b);
Y = [];
fac_size = zeros(options.AL_maxiter, 1);
diag_index = (1:(n+1):n^2).';
previous_pinf = inf;
grad_floor = min(options.tolgradnorm, 0.1*options.tol);
inner_iter = options.TR_maxiter;
inner_calls = options.TR_maxinner;

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
% The user budgets are starting values. Spend extra inner work only if an
% inaccurate subproblem also stalls feasibility, then restore the budgets.

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    fac_size(iter) = p;
    problem.M = obliquefactoryNTrans(p, n);
    dual_c = c - At*y;
    opts.tolgradnorm = max(grad_floor, min(tau, 0.1*previous_pinf));
    opts.maxiter = inner_iter;
    opts.maxinner = min(inner_calls, max(1, problem.M.dim()));
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    x = X(:);
    obj = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;     
    y = y - sigma*Axb;
    eS = reshape(c - At*y, n, n);
    z = sum(X.*eS).';
    S = eS;
    S(diag_index) = S(diag_index) - z;
    [vS, dS] = eig(S, 'vector');
    [dS, order] = sort(dS);
    vS = vS(:, order);
    dinf = max(0, -dS(1))/(1+max(0, dS(end)));
    by = b'*y + sum(z);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    eta = max([gap, pinf, dinf]);
    % A large theta must not discard significant modes. Tighten the cutoff
    % further near convergence to avoid repeated loss of feasible digits.
    rank_tol = min([options.theta, 1e-3, ...
                   0.1*sqrt(max(options.tol, min(1, eta)))]);
    r = max(1, sum(e >= rank_tol*e(1)));
    if options.verbosity
        fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,       toc(timespend));
    end
    if eta < options.tol
        if options.verbosity, fprintf('Optimality is reached!\n'); end
        break;
    end
    if ~isfinite(eta) || ~isfinite(gradnorm)
        data.status = 2;
        break;
    end
    inner_inaccurate = gradnorm > max(opts.tolgradnorm, 0.1*norm(Axb));
    if inner_inaccurate && pinf > 0.9*previous_pinf ...
            && info(end).iter >= opts.maxiter
        inner_iter = min(max(8, options.TR_maxiter), 2*inner_iter);
        inner_calls = min(max(30, options.TR_maxinner), 2*inner_calls);
    elseif gradnorm <= opts.tolgradnorm || pinf < 0.7*previous_pinf
        inner_iter = max(options.TR_maxiter, ceil(inner_iter/2));
        inner_calls = max(options.TR_maxinner, ceil(inner_calls/2));
    end
    if r <= p - 1
        Y = V(:,1:r)'.*e(1:r);
        Y = Y./sqrt(sum(Y.^2, 1));
        p = r;
    end
    tau = max(0.4*tau, grad_floor);
    nne = min([sum(dS < -tau), options.delta, n-p]);
    if nne > 0
        U = [zeros(p, n); vS(:,1:nne)'];
        Y = [Y; zeros(nne,n)]; %#ok<AGROW> Staircase increases the factor size.
        % The dual certificate uses y - sigma*(A(X)-b), so its slack is
        % precisely the current augmented objective's gradient matrix.
        % dual_c retains the just-solved subproblem's multipliers: search
        % that same objective before moving on to the next subproblem.
        Y = line_search(Y, U);
        p = p + nne;
    end
    % A dead band prevents alternating penalty changes on nearly equal
    % residuals; extra inner work is reserved for measured stagnation.
    if pinf < options.kappa1*gradnorm
        sigma = max(sigma/gamma, options.sigma_min);
    elseif pinf > options.kappa2*gradnorm
        sigma = min(sigma*gamma, options.sigma_max);
    end
    previous_pinf = pinf;
end
data.X = X;
data.y = y;
data.S = S;
data.z = z;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
data.fac_size = fac_size(1:iter);
if data.status == 0 && eta >= options.tol
    data.status = 1;
    if options.verbosity, fprintf('Iteration maximum is reached!\n'); end
end

if options.verbosity
    fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
end

    function nY = line_search(Y, U)
        % The staircase direction satisfies Y'*U = 0. The normalized
        % pullback is even in alpha; its derivative in theta = alpha^2
        % is half the Riemannian Hessian curvature along U.
        nY = Y;
        ls_ynorm = sqrt(sum(Y.^2, 1));
        if any(~isfinite(ls_ynorm)) || any(ls_ynorm == 0)
            return;
        end
        % Restore unit columns after possible rank truncation.
        Y = Y./ls_ynorm;
        nY = Y;
        ls_u2 = sum(U.^2, 1);
        if any(~isfinite(ls_u2)) || max(ls_u2) == 0
            return;
        end

        ls_X0 = Y'*Y;
        ls_diag_index = 1:(size(ls_X0, 1)+1):numel(ls_X0);
        ls_X0(ls_diag_index) = 1;
        ls_UtU = U'*U;
        ls_dX = ls_UtU ...
            - 0.5*(ls_u2.'.*ls_X0 + ls_X0.*ls_u2);
        ls_dX(ls_diag_index) = 0;
        ls_residual = A*ls_X0(:) - b;
        ls_linear = dual_c + sigma*(A'*ls_residual);
        ls_slope = ls_linear'*ls_dX(:);
        if ~isfinite(ls_slope) || ls_slope >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(max(ls_u2)));
        if options.line_search == 0
            ls_alpha = min(options.alpha, 1/sqrt(max(ls_u2)));
        end
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            ls_scale = sqrt(1 + ls_theta*ls_u2);
            % Stable evaluation of sqrt(1 + theta*u2) - 1.
            ls_shift = (ls_theta*ls_u2)./(1 + ls_scale);
            % Exact change in the column-normalized Gram matrix.
            ls_deltaX = (ls_theta*ls_UtU ...
                - (ls_shift.' + ls_shift + ls_shift.'*ls_shift).*ls_X0) ...
                ./ (ls_scale.'*ls_scale);
            ls_deltaX(ls_diag_index) = 0;
            ls_AdeltaX = A*ls_deltaX(:);
            ls_change = ls_linear'*ls_deltaX(:) ...
                + 0.5*sigma*(ls_AdeltaX'*ls_AdeltaX);
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                nY = Y + ls_alpha*U;
                nY = nY./sqrt(sum(nY.^2, 1));
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep the normalized base point if no trial is acceptable.
    end

    function [f, store] = cost(Y, store)
        gram = Y'*Y;
        store.residual = A*gram(:) - b;
        % Drop the constant ||y||^2/(2*sigma) to avoid cancellation when
        % comparing trust-region costs near convergence.
        f = dual_c'*gram(:) + 0.5*sigma*(store.residual'*store.residual);
    end

    function [G, store] = grad(Y, store)
        if ~isfield(store, 'residual')
            [~, store] = cost(Y, store);
        end
        store.slack = reshape(dual_c + sigma*(At*store.residual), n, n);
        eG = 2*Y*store.slack;
        store.YeG = sum(Y.*eG, 1);
        G = eG - Y.*store.YeG;
    end

    function [H, store] = hess(Y, U, store)
        YU = Y'*U;
        AyU = reshape(A'*(At'*YU(:)), n, n);
        if ~isfield(store, 'slack')
            [~, store] = grad(Y, store);
        end
        eH = 2*U*store.slack + 4*sigma*(Y*AyU);
        H = eH - Y.*sum(Y.*eH, 1) - U.*store.YeG;
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U, 1);
        M.tangent = @(X, U) U - X.*sum(X.*U, 1);

        M.retr = @retraction;
        function y = retraction(x, d)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2, 1));
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d, 1);

        function x = random(n, m)
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end
