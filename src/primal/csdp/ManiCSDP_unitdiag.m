% This function solves CSDPs with unit-diagonal using the primal approach:
% Min  <C, X>
% s.t. A(X) = b
%      X >= 0
%      diag(X) = 1

function [X, obj, data] = ManiCSDP_unitdiag(At, b, c, K, options)
% Defaults selected on example_bqp_complex with d=20,30 and tol=1e-8:
% sigma0=2e-2, sigma_min=1e-2, gamma=2, p0=2, theta=1e-3,
% delta=8, TR_maxiter=4, TR_maxinner=25, kappa1=kappa2=1, line_search=1.
% Inner budgets adapt on stagnation. Both staircase modes enforce decrease;
% line_search=0 uses alpha as its initial trial step.

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s) ...
        || ~isreal(K.s) || ~isfinite(K.s) || K.s < 1 || K.s ~= fix(K.s)
    error('K.s must be a positive scalar integer.');
end

n = K.s;
if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gamma')
    if isfield(options,'gama'), options.gamma = options.gama;
    else, options.gamma = 2; end
end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 2e-2; end
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

if options.verbosity
    fprintf('ManiCSDP is starting...\n');
    fprintf('CSDP size: n = %i, m = %i\n', n, numel(b));
end
warning_state = warning('off', 'manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));

% Only Hermitian parts contribute on Hermitian X. Canonicalize once so
% every Hessian application can use the operator without dense symmetrizing.
transpose_index = reshape(1:n*n, n, n).';
At = 0.5*(At + conj(At(transpose_index(:), :)));
C = reshape(c, n, n);
C = (C+C')/2;
c = C(:);
A = At';
p = min(n, max(1, round(options.p0)));
% Preserve the legacy convention allowing sigma0 below sigma_min;
% sigma_min bounds subsequent penalty decreases.
sigma = min(options.sigma_max, max(realmin, options.sigma0));
gamma = options.gamma;
tau = options.tau;
y = zeros(length(b), 1);
normb = 1 + norm(b);
Y    = [];
fac_size = zeros(options.AL_maxiter, 1);
diag_index = (1:n+1:n*n).';
previous_pinf = inf;
grad_floor = min(options.tolgradnorm, 0.1*options.tol);
inner_iter = options.TR_maxiter;
inner_calls = options.TR_maxinner;

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    fac_size(iter) = p;
    problem.M = obliquefactorycomplexNTrans(p, n);
    dual_c = c - At*y;
    opts.tolgradnorm = max(grad_floor, min(tau, 0.1*previous_pinf));
    opts.maxiter = inner_iter;
    opts.maxinner = min(inner_calls, problem.M.dim());
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    x = X(:);
    obj = real(c'*x);
    Axb = real(A*x - b);
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    eS = reshape(c - At*y, n, n);
    eS = (eS+eS')/2;
    z = real(sum(conj(X).*eS, 1));
    S = eS;
    S(diag_index) = S(diag_index) - z.';
    if ~isfinite(obj) || ~isfinite(pinf) || ~isfinite(gradnorm) ...
            || any(~isfinite(S(:)))
        gap = inf; dinf = inf; eta = inf;
        data.status = 2;
        break;
    end
    [vS, dS] = eig(S, 'vector');
    [dS, order] = sort(real(dS));
    vS = vS(:, order);
    dinf = max(0, -dS(1))/(1+max(0, dS(end)));
    by = real(b'*y) + sum(z);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    eta = max([gap, pinf, dinf]);
    % Cap aggressive truncation and tighten it near the requested accuracy.
    rank_tol = min([options.theta, 1e-3, sqrt(max(options.tol, min(1, eta)))]);
    r = max(1, sum(e >= rank_tol*e(1)));
    if options.verbosity
        fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,       toc(timespend));
    end
    if eta < options.tol
        if options.verbosity, fprintf('Optimality is reached!\n'); end
        break;
    end
    % Refine only when both feasibility and the inner solve stall, rather
    % than spending a large fixed trust-region budget on every outer pass.
    inner_inaccurate = gradnorm > max(opts.tolgradnorm, 0.1*norm(Axb));
    if inner_inaccurate && pinf > 0.9*previous_pinf ...
            && info(end).iter >= opts.maxiter
        inner_iter = min(max(4, options.TR_maxiter), 2*inner_iter);
        inner_calls = min(max(30, options.TR_maxinner), 2*inner_calls);
    elseif gradnorm <= opts.tolgradnorm || pinf < 0.7*previous_pinf
        inner_iter = max(options.TR_maxiter, ceil(inner_iter/2));
        inner_calls = max(options.TR_maxinner, ceil(inner_calls/2));
    end
    if r <= p - 1
        Y = V(:,1:r)'.*e(1:r);
        Y = Y./sqrt(real(sum(conj(Y).*Y, 1)));
        p = r;
    end
    % Complex rank-one factors still have phase directions. Add rows only
    % for genuine negative curvature; never exceed the matrix dimension.
    tau = max(0.4*tau, grad_floor);
    nne = min([sum(dS < -tau), options.delta, n-p]);
    if nne > 0
        U = [zeros(p,n); vS(:,1:nne)'];
        Y = [Y; zeros(nne,n)]; %#ok<AGROW> The staircase increases the rank.
        % S is the slack of the just-solved augmented objective. Search
        % that same objective (dual_c), before starting the next subproblem.
        Y = line_search(Y, U);
        p = p + nne;
    end
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
    fprintf('ManiCSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
end

    function nY = line_search(Y, U)
        % The complex staircase direction satisfies Y'*U = 0. The normalized
        % pullback is even in alpha; its derivative in theta = alpha^2
        % is half the Riemannian Hessian curvature along U.
        nY = Y;
        ls_ynorm = sqrt(real(sum(conj(Y).*Y, 1)));
        if any(~isfinite(ls_ynorm)) || any(ls_ynorm == 0)
            return;
        end
        % Restore unit columns after possible rank truncation.
        Y = Y./ls_ynorm;
        nY = Y;
        ls_u2 = real(sum(conj(U).*U, 1));
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
        ls_residual = real(A*ls_X0(:)) - b;
        ls_linear = dual_c + sigma*(At*ls_residual);
        ls_slope = real(ls_linear'*ls_dX(:));
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
            ls_AdeltaX = real(A*ls_deltaX(:));
            ls_change = real(ls_linear'*ls_deltaX(:)) ...
                + 0.5*sigma*real(ls_AdeltaX'*ls_AdeltaX);
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                nY = Y + ls_alpha*U;
                nY = nY./sqrt(real(sum(conj(nY).*nY, 1)));
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep the normalized base point if no trial is acceptable.
    end

    function [f, store] = cost(Y, store)
        gram = Y'*Y;
        store.residual = real(A*gram(:)) - b;
        % Omit the constant ||y||^2/(2*sigma) to keep cost differences
        % accurate near convergence. Cache each point independently.
        f = real(dual_c'*gram(:)) + 0.5*sigma*(store.residual'*store.residual);
    end

    function [G, store] = grad(Y, store)
        if ~isfield(store, 'residual')
            [~, store] = cost(Y, store);
        end
        store.eS = reshape(dual_c + sigma*(At*store.residual), n, n);
        eG = 2*Y*store.eS;
        store.YeG = real(sum(conj(Y).*eG, 1));
        G = eG - Y.*store.YeG;
    end

    function [H, store] = hess(Y, U, store)
        YU  = Y'*U;
        AyU = reshape(At*(real(At'*YU(:))), n, n);
        if ~isfield(store, 'eS') || ~isfield(store, 'YeG')
            [~, store] = grad(Y, store);
        end
        eH  = 2*U*store.eS + 4*sigma*(Y*AyU);
        H   = eH - Y.*real(sum(conj(Y).*eH,1)) - U.*store.YeG;
    end

    function M = obliquefactorycomplexNTrans(n, m)
        M.dim = @() (2*n-1)*m;
        M.inner = @(x, d1, d2) real(d1(:)'*d2(:));
        M.norm = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*real(sum(conj(X).*U,1));
        M.tangent = @(X, U) U - X.*real(sum(conj(X).*U,1));

        M.retr = @retraction;
        function y = retraction(x, d)
            xtd = x + d;
            y = xtd./sqrt(real(sum(conj(xtd).*xtd,1)));
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*real(sum(conj(x2).*d,1));

        function x = random(n, m)
            x = randn(n, m)+1i*randn(n,m);
            x = x./sqrt(real(sum(conj(x).*x, 1)));
        end
    end
end
