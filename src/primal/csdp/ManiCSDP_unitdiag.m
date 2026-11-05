% This function solves CSDPs with unit-diagonal using the primal approach:
% Min  <C, X>
% s.t. A(X) = b
%      X >= 0
%      diag(X) = 1

function [X, obj, data] = ManiCSDP_unitdiag(At, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s)
    error('K.s must be a positive scalar integer.');
end

n = K.s;
if ~isfield(options,'p0'); options.p0 = 2; end
% if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-3; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiCSDP is starting...\n');
fprintf('CSDP size: n = %i, m = %i\n', n, size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');

A = At';
p = options.p0;
sigma = options.sigma0;
gamma = options.gamma;
tau = options.tau;
y = zeros(length(b), 1);
normb = 1 + norm(b);
Y    = [];
U    = [];
Axb  = [];
fac_size = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    fac_size = [fac_size; p];
    problem.M = obliquefactorycomplexNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    x = X(:);
    obj = real(c'*x);
    Axb = real(A*x - b);
    pinf = norm(Axb)/normb;
    y = y - sigma*real(Axb);
    eS = reshape(c - At*y, n, n);
    eS = (eS+eS')/2;
    z = real(sum(conj(X).*eS, 1));
    S = eS - diag(z);
    [vS, dS] = eig(S, 'vector');
    dS = sort(real(dS));
    dinf = max(0, -dS(1))/(1+abs(dS(end)));
    by = real(b'*y) + sum(z);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    r = sum(e >= options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,       toc(timespend));
    eta = max([gap, pinf, dinf]);
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 20) == 0
        if iter > 50 && gap > gap0 && pinf > pinf0 && dinf > dinf0
            data.status = 2;
            fprintf('Slow progress!\n');
            break;
        else
            gap0 = gap;
            pinf0 = pinf;
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
       Y = Y./sqrt(real(sum(conj(Y).*Y,1)));
    end
    if pinf < options.kappa1*gradnorm
          sigma = max(sigma/gamma, options.sigma_min);
    elseif pinf > options.kappa2*gradnorm
          sigma = min(sigma*gamma, options.sigma_max);
    end
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
data.fac_size = fac_size;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiCSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

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
        ls_residual = A*ls_X0(:) - b - y/sigma;
        ls_linear = c + sigma*(A'*ls_residual);
        ls_slope = real(ls_linear'*ls_dX(:));
        if ~isfinite(ls_slope) || ls_slope >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(max(ls_u2)));
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
        X   = Y'*Y;
        x   = X(:);
        Axb = A*x - b - y/sigma;
        f   = real(c'*x + 0.5*sigma*(Axb'*Axb));
        store.eS = reshape(c + sigma*At*Axb, n, n);
        store.eS = (store.eS + store.eS') / 2;
    end

    function [G, store] = grad(Y, store)
        eG = 2*Y*store.eS;
        store.YeG = real(sum(conj(Y).*eG, 1));
        G = eG - Y.*store.YeG;
    end

    function [H, store] = hess(Y, U, store)
        YU  = Y'*U;
        AyU = reshape(At*(real(At'*YU(:))), n, n);
        eH  = 2*U*store.eS + 4*sigma*(Y*AyU);
        H   = eH - Y.*real(sum(conj(Y).*eH,1)) - U.*store.YeG;
    end

    function M = obliquefactorycomplexNTrans(n, m)
        M.dim = @() (2*n-1)*m;
        M.inner = @(x, d1, d2) real(sum(conj(d1).*d2,'all'));
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
