% This function solves SDPs using the dual approach:
%    sup  <C, X> + <c, w>
%    s.t. A(X) + B(w) = b
%         X in S_+^{n}
%         w in R^l

function [X, obj, data] = ManiDSDP(A, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s)
    error('K.s must be a positive scalar integer.');
end
if ~isfield(K, 'f') || isempty(K.f)
    K.f = 0;
end

% if ~isfield(options,'p0'); options.p0 = 1; end
if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter = 1000; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.01; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 0.1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1; end
if ~isfield(options,'line_search'); options.line_search = 1; end

n = K.s;
fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));
warning_state = warning('off', 'manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));

normc = 1 + norm(c);
B = A(:, 1:K.f);
A = A(:, K.f+1:end);
cf = c(1:K.f);
c = c(K.f+1:end);

if ~isfield(options,'dAAt')
    iA = ((A*A')\A)'; % (inv(A*A')*A)'
else
    iA = (sparse(1:size(b,1),1:size(b,1),options.dAAt)\A)';
end
bA = iA*b;
iAB = iA*B;

p = options.p0;
sigma = options.sigma0; 
gamma = options.gamma;
tau = options.tau;
x = zeros(n^2, 1);
w = zeros(K.f, 1);
Y = [];
U = [];

problem.costgrad = @costgrad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = options.TR_maxinner;    % maximum Hessian calls per iteration
opts.maxiter   = options.TR_maxiter;
% opts.tolgradnorm = options.tolgradnorm;
% opts.trscache = false;
data.status = 0;
timespend = tic;
for iter = 1:options.ADMM_maxiter
    problem.M = euclideanfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    S = Y*Y';
    sc = S(:) - c;
    y = iA'*sc;
    As = A'*y - sc;
    Af = B'*y - cf;
    pinf = (norm(As)+norm(Af))/normc;
    by = b'*y;
    x = x + sigma*(iAB*(Af - w/sigma) + A'*(iA'*(As - x/sigma)) - As);
    w = w - sigma*Af;
    X = reshape(x + bA, n, n);
    [vX, dX] = eig(X, 'vector');
    obj = c'*(x + bA) + cf'*w;
    dinf = max(0, -dX(1))/(1+abs(dX(end)));
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    [V, D, ~] = svd(Y, 'econ');
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma, toc(timespend));
    eta = max([gap, pinf, dinf]);
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 50) == 0
        if iter > 100 && gap > gap0 && pinf > pinf0 && dinf > dinf0
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
    nne = min(sum(dX < -tau), options.delta);
    if options.line_search == 1
        U = [zeros(n, p) vX(:,1:nne)];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y zeros(n, nne)];
    else
        Y = [Y options.alpha*vX(:,1:nne)];
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
data.w = w;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
        
    function val = co(Y)
        S = Y*Y';
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        val = b'*y + 0.5*sigma*(As'*As + Af'*Af);
    end

    function nY = line_search(Y, U)
        % The staircase direction satisfies Y*U' = 0. With theta =
        % alpha^2, S(theta) = Y*Y' + theta*U*U', so the derivative in
        % theta is half the Hessian curvature along U.
        nY = Y;
        ls_u2 = sum(U.^2, 'all');
        if ~isfinite(ls_u2) || ls_u2 == 0 || any(~isfinite(Y), 'all')
            return;
        end

        % Use local names to preserve the Manopt callbacks' shared state.
        ls_S0 = Y*Y';
        ls_sc0 = ls_S0(:) - c;
        ls_y0 = iA'*ls_sc0;
        ls_As0 = A'*ls_y0 - ls_sc0 - x/sigma;
        ls_Af0 = B'*ls_y0 - cf - w/sigma;
        ls_dS = U*U';
        ls_dy = iA'*ls_dS(:);
        ls_dAs = A'*ls_dy - ls_dS(:);
        ls_dAf = B'*ls_dy;
        ls_slope = b'*ls_dy ...
            + sigma*(ls_As0'*ls_dAs + ls_Af0'*ls_dAf);
        ls_quadratic = 0.5*sigma*(ls_dAs'*ls_dAs + ls_dAf'*ls_dAf);
        if ~isfinite(ls_slope) || ls_slope >= 0 ...
                || ~isfinite(ls_quadratic)
            return;
        end

        ls_alpha = min(1, 1/sqrt(ls_u2));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            % Exact dual-objective change, without subtracting close costs.
            ls_change = ls_theta*(ls_slope + ls_theta*ls_quadratic);
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                nY = Y + ls_alpha*U;
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep Y if no trial satisfies second-order Armijo decrease.
    end

    function [f, G, store] = costgrad(Y, store)
        S = Y*Y';
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        f = b'*y + 0.5*sigma*(As'*As + Af'*Af);
        X = reshape(bA + sigma*(iAB*Af + A'*(iA'*As) - As), n, n);
        G = 2*X*Y;
    end

    function [H, store] = hess(Y, U, store)      
        YU = U*Y';
        yAU = reshape((YU(:)'*iA)*A, n, n);
        H = 2*X*U + 2*sigma*(U*(Y'*Y) + Y*(U'*Y)) + 4*sigma*(reshape((YU(:)'*iAB)*iAB' + (yAU(:)'*iA)*A, n, n) - 2*yAU)*Y;
    end
end