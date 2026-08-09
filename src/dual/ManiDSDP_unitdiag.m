% This function solves SDPs with dual unit-diagonal (diag(S) = 1) using the dual approach:
%    sup  <C, X> + <c, w>
%    s.t. A(X) + B(w) = b
%         X in S_+^{n}
%         w in R^l
%         diag(S) = 1

function [X, obj, data] = ManiDSDP_unitdiag(A, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(K) || ~isfield(K, 's') || ~isscalar(K.s)
    error('K.s must be a positive scalar integer.');
end
if ~isfield(K, 'f') || isempty(K.f)
    K.f = 0;
end

% if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter = 300; end
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
if ~isfield(options,'kappa1'); options.kappa1 = 1e1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1e2; end
if ~isfield(options,'line_search'); options.line_search = 0; end

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

if ~isfield(options, 'dAAt')
    iA = ((A*A')\A)'; % (inv(A*A')*A)'
else
    iA = (sparse(1:size(b,1),1:size(b,1),options.dAAt)\A)';
end
bA = iA*b;
diag_index = (1:(n+1):n^2).';

p = options.p0;
sigma = options.sigma0;
gamma = options.gamma;
tau = options.tau;
x = zeros(n^2, 1);
w = zeros(K.f, 1);
eG = zeros(p, n);
Y = [];
U = [];
fac_size = [];
seta = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = options.TR_maxinner;    % maximum Hessian calls per iteration
opts.maxiter   = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;
% opts.trscache = false;
data.status = 0;
timespend = tic;
for iter = 1:options.ADMM_maxiter
    fac_size = [fac_size; p];
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    S = Y'*Y;
    sc = S(:) - c;
    y = iA'*sc;
    As = A'*y - sc;
    Af = B'*y - cf;
    pinf = (norm(As)+norm(Af))/normc;
    by = b'*y;
    x = x - sigma*As;
    w = w - sigma*Af;
    eX = reshape(x + bA, n, n);
    z = sum(S .* eX, 1).';
    X = eX;
    X(diag_index) = X(diag_index) - z;
    [vX, dX] = eig(X, 'vector');
    % [vX, dX] = eigs(X, options.delta, 'smallestreal');
    % dX = diag(dX);
    obj = c'*eX(:) + cf'*w + sum(z);
    dinf = max(0, -dX(1))/(1+abs(dX(end)));
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,   toc(timespend));
    eta = max([gap, pinf, dinf]);
    seta = [seta; eta];
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
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    if p == 1
       nne = 1;
    else
        tau = max(0.4*tau, 1e-8);
        nne = min(sum(dX < -tau), options.delta);
    end
    if options.line_search == 1
        U = [zeros(p, n); vX(:,1:nne)'];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y; zeros(nne, n)];
    else
        Y = [Y; options.alpha*vX(:,1:nne)'];
        Y = Y./sqrt(sum(Y.^2));
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
data.fac_size = fac_size;
data.seta = seta;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
       

    function nY = line_search(Y, U)
        % Negative-curvature line search on the oblique manifold.
        % For the staircase direction U, Y'*U = 0, so the retracted
        % pullback is even in alpha. With theta = alpha^2, its derivative
        % at theta = 0 is one half of the Riemannian Hessian curvature.
        % Hence, use a second-order Armijo test in theta rather than a
        % first-order Armijo test in alpha.

        nY = Y;
        armijo = 0.1;
        contraction = 0.5;       % Backtracking factor for theta.
        max_backtracks = 25;

        % Keep the base point on the oblique manifold (unit-norm columns
        % in this transposed representation).
        ynorm = sqrt(sum(Y.^2, 1));
        if any(~isfinite(ynorm)) || any(ynorm <= sqrt(eps))
            return;
        end
        Y = bsxfun(@rdivide, Y, ynorm);
        nY = Y;

        u2 = sum(U.^2, 1);
        max_u2 = max(u2);
        if ~isfinite(max_u2) || max_u2 <= eps
            return;
        end

        % The staircase construction requires the stronger orthogonality
        % Y'*U = 0, not merely columnwise tangency.
        orth_scale = max(1, norm(Y, 'fro')*norm(U, 'fro'));
        if norm(Y'*U, 'fro') > 100*eps*orth_scale
            return;
        end

        % Base objective data. Unique local names avoid changing the
        % shared As/Af variables used by the Manopt callbacks.
        ls_S0 = Y'*Y;
        ls_sc0 = ls_S0(:) - c;
        ls_y0 = iA'*ls_sc0;
        ls_As0 = A'*ls_y0 - ls_sc0 - x/sigma;
        ls_Af0 = B'*ls_y0 - cf - w/sigma;
        linear0 = b'*ls_y0;
        quad0 = 0.5*sigma*(ls_As0'*ls_As0 + ls_Af0'*ls_Af0);
        cost0 = linear0 + quad0;

        % Exact derivative of phi(theta) = co(Retr_Y(sqrt(theta)U)) at
        % theta = 0. Since Y'*U = 0,
        %   S'(0) = U'*U - 1/2*(Diag(u2)S0 + S0*Diag(u2)).
        ls_Sdot = U'*U ...
            - 0.5*(bsxfun(@times, u2.', ls_S0) ...
                 + bsxfun(@times, ls_S0, u2));
        ls_dsc = ls_Sdot(:);
        ls_dy = iA'*ls_dsc;
        ls_dAs = A'*ls_dy - ls_dsc;
        ls_dAf = B'*ls_dy;
        slope = b'*ls_dy ...
              + sigma*(ls_As0'*ls_dAs + ls_Af0'*ls_dAf);

        % A proper negative-curvature step requires phi'(0) < 0, or
        % equivalently <U, Hess co(Y)[U]> = 2*phi'(0) < 0.
        scale0 = max([1, abs(linear0), abs(quad0), abs(cost0)]);
        curvature_tol = 100*eps*scale0;
        if ~isfinite(slope) || slope >= -curvature_tol
            return;
        end

        % Start with alpha <= 1 and cap every column's initial rotation at
        % 45 degrees: alpha^2*||U(:,j)||^2 <= 1.
        theta = min(1, 1/max_u2);

        for j = 1:max_backtracks
            predicted_decrease = -theta*slope;
            if predicted_decrease <= 1000*eps*scale0
                return;
            end

            alpha = sqrt(theta);
            trialY = Y + alpha*U;
            trialY = bsxfun(@rdivide, trialY, sqrt(sum(trialY.^2, 1)));

            % Compute the objective difference from the exact quadratic
            % expansion in S, avoiding cancellation in co(trialY)-co(Y).
            ls_dsc_trial = trialY'*trialY - ls_S0;
            ls_dsc_trial = ls_dsc_trial(:);
            ls_dy_trial = iA'*ls_dsc_trial;
            ls_dAs_trial = A'*ls_dy_trial - ls_dsc_trial;
            ls_dAf_trial = B'*ls_dy_trial;
            delta_cost = b'*ls_dy_trial ...
                + sigma*(ls_As0'*ls_dAs_trial + ls_Af0'*ls_dAf_trial) ...
                + 0.5*sigma*(ls_dAs_trial'*ls_dAs_trial ...
                           + ls_dAf_trial'*ls_dAf_trial);

            % Second-order Armijo condition:
            %   phi(theta) <= phi(0) + armijo*theta*phi'(0).
            rhs = armijo*theta*slope;
            roundoff = 10*eps*max([1, abs(cost0), abs(delta_cost)]);
            if isfinite(delta_cost) && delta_cost < 0 ...
                    && delta_cost <= rhs + roundoff
                nY = trialY;
                return;
            end

            theta = contraction*theta;
        end
    end

      function [f, store] = cost(Y, store)
        S = Y'*Y;
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        f = b'*y + 0.5*sigma*(As'*As + Af'*Af);
      end

      function [G, store] = grad(Y, store)
        X = reshape(bA - sigma*As, n, n);
        eG = 2*Y*X;
        G = eG - Y.*sum(Y.*eG);
      end

    function [H, store] = hess(Y, U, store)
        YU = Y'*U;
        yAU = reshape(A'*(iA'*YU(:)), n, n);
        eH = 2*U*X - 4*sigma*(Y*yAU) + 2*sigma*((Y*U')*Y+(Y*Y')*U);
        H = eH - Y.*sum(Y.*eH) - U.*sum(Y.*eG);
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
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
            x = x./sqrt(sum(x.^2));
        end
    end
end

