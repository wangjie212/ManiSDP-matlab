% This function solves CSDPs with unit trace using the primal approach:
% Min  <C, X>
% s.t. A(X) = b
%      X >= 0
%      tr(X) = 1

function [X, obj, data] = ManiCSDP_unittrace(At, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s)
    error('K.s must be a positive scalar integer.');
end

n = K.s;
if ~isfield(options,'p0'); options.p0 = 1; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 1000; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.05; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 40; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 3; end
if ~isfield(options,'kappa1'); options.kappa1 = 1e-5; end
if ~isfield(options,'kappa2'); options.kappa2 = 1e-4; end
if ~isfield(options,'line_search'); options.line_search = 1; end

fprintf('ManiCSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));

A = At';
p = options.p0;
sigma = options.sigma0;
gamma = options.gamma;
tau = options.tau;
y = zeros(length(b),1);
normb = 1 + norm(b);
if isfield(options, 'Y0') 
    Y = options.Y0; 
else
    Y = [];
end
U = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;
opts.maxinner = options.TR_maxinner;
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    problem.M = spherecomplexfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y*Y';
    x = X(:);
    obj = real(c'*x);
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    eS = reshape(c - At*y, n, n);
    eS = (eS + eS')/2;
    z = real(sum(conj(eS).*X, 'all'));
    S = eS - z*eye(n);
    [vS, dS] = eig(S, 'vector');
    dinf = max(0, -dS(1))/(1+dS(end));
    by = real(b'*y) + z;
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [V, D, ~] = svd(Y, 'econ');
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e >= options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,    r,    p,    sigma,   toc(timespend));
    eta = max([pinf, gap, dinf]);
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
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    tau = max(0.4*tau, 1e-8);
    nne = min(sum(dS < -tau), options.delta);
    if options.line_search == 1
       U = [zeros(n, p) vS(:,1:nne)];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y zeros(n, nne)];
    else
        Y = [Y options.alpha*vS(:,1:nne)];
        Y = Y/norm(Y, 'fro');
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
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiCSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    function val = co(Y)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = real(c'*x + sigma/2*(Axb'*Axb));
    end

    function nY = line_search(Y, U)
        % The complex staircase direction satisfies Y*U' = 0, so the pullback
        % along the sphere retraction is even in alpha. Its derivative
        % in theta = alpha^2 is half the Riemannian Hessian curvature.
        nY = Y;
        ls_ynorm = norm(Y, 'fro');
        if ~isfinite(ls_ynorm) || ls_ynorm == 0
            return;
        end
        % Restore unit trace after possible rank truncation.
        Y = Y/ls_ynorm;
        nY = Y;
        ls_u2 = real(sum(conj(U).*U, 'all'));
        if ~isfinite(ls_u2) || ls_u2 == 0
            return;
        end

        ls_X0 = Y*Y';
        ls_dX = U*U' - ls_u2*ls_X0;
        ls_residual = A*ls_X0(:) - b - y/sigma;
        ls_AdX = A*ls_dX(:);
        ls_slope = real(c'*ls_dX(:) + sigma*(ls_residual'*ls_AdX));
        ls_quadratic = 0.5*sigma*real(ls_AdX'*ls_AdX);
        if ~isfinite(ls_slope) || ls_slope >= 0 ...
                || ~isfinite(ls_quadratic)
            return;
        end

        ls_alpha = min(1, 1/sqrt(ls_u2));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            ls_t = ls_theta/(1 + ls_theta*ls_u2);
            % Retracted X = X0 + t*dX. This gives the exact change in
            % the augmented Lagrangian without subtracting close costs.
            ls_change = ls_t*ls_slope + ls_t^2*ls_quadratic;
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                nY = Y + ls_alpha*U;
                nY = nY/norm(nY, 'fro');
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep the normalized base point if no trial is acceptable.
    end

    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = real(c'*x + sigma/2*(Axb'*Axb));
        eS = reshape(c+sigma*At*Axb, n, n);
        eS=(eS+eS')/2;
        store.z = real(sum(conj(X).*eS, 'all'));
        store.G = 2*eS*Y - 2*store.z*Y;
        store.eS = eS;
    end
    
    function [G, store] = grad(~, store)
        G = store.G;
    end

    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU= reshape(YU(:)'*At*A, n, n);
        H = 2*store.eS*U + 4*sigma*(AyU*Y);
%        H = H - sum(H.*Y,'all')*Y - 2*store.z*U;
        H = H - real(trace(H*Y'))*Y - 2*store.z*U;
    end
end

