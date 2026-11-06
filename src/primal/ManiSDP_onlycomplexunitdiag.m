% Solves complex unit-diagonal SDP via paired-oblique manifold on real conversion.
%
% Problem (real, after convertCtoR):
%   min  <C, X>
%   s.t. X >= 0
%        X(i,i) + X(n+i,n+i) = 2,  i = 1,...,n
%
% Factorization: X = Y'*Y, where Y is (p x 2n).
% Constraint encoded implicitly: ||col_i(Y)||^2 + ||col_{n+i}(Y)||^2 = 2.
%
% Analogous to ManiSDP_onlyunitdiag for real MaxCut; no explicit At matrix needed.
% Input C must be the 2n x 2n real matrix from convertCtoR (already divided by 2).

function [X, obj, data] = ManiSDP_onlycomplexunitdiag(C, options)

if nargin < 2 || isempty(options)
    options = struct();
end

n2 = size(C, 1);   % = 2n
n  = n2 / 2;

if ~isfield(options,'p0');          options.p0          = 2;    end
if ~isfield(options,'AL_maxiter');  options.AL_maxiter  = 20;   end
if ~isfield(options,'tol');         options.tol         = 1e-8; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'theta');       options.theta       = 1e-1; end
if ~isfield(options,'delta');       options.delta       = 8;    end
if ~isfield(options,'alpha');       options.alpha       = 0.5;  end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 100;  end
if ~isfield(options,'TR_maxiter');  options.TR_maxiter  = 40;   end
if ~isfield(options,'line_search'); options.line_search = 0;    end

fprintf('ManiSDP_onlycomplexunitdiag is starting...\n');
fprintf('Complex SDP size: n = %d, real size: 2n = %d\n', n, n2);

p       = options.p0;
tau = options.tau;
Y       = [];
U       = [];
YC      = [];      % shared closure: set in cost(), used in grad()
eG_pair = [];      % shared closure: set in cost(), used in grad() and hess()

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity   = 0;
opts.maxinner    = options.TR_maxinner;
opts.maxiter     = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    problem.M = pairedobliquefactory(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;

    % Compute objective and dual certificate
    X        = Y'*Y;
    Y1       = Y(:, 1:n);
    Y2       = Y(:, n+1:n2);
    YC_out   = Y * C;
    eG1      = sum(Y1 .* YC_out(:, 1:n));     % 1 x n
    eG2      = sum(Y2 .* YC_out(:, n+1:n2));  % 1 x n
    eG_pair  = (eG1 + eG2) / 2;               % 1 x n, dual multipliers y_i
    obj      = 2 * sum(eG_pair);              % = trace(C * X)

    % Dual slack: S = C - diag([y; y]), where y_i = 2*eG_pair_i
    % (from KKT: (Y*C)[:,i] = y_i/2 * Y[:,i], so y_i = 2*eG_pair_i,
    %  and S = C - sum_i y_i*(e_i*e_i'+e_{n+i}*e_{n+i}')/2
    %          = C - diag([eG_pair, eG_pair]))
    z  = [eG_pair, eG_pair];          % 1 x 2n
    S  = C - diag(z(:));            % 2n x 2n
    [vS, dS] = eig(full(S), 'vector');
    dS = sort(real(dS));
    dinf = max(0, -dS(1)) / (1 + abs(dS(end)));

    [~, D, V] = svd(Y, 'econ');
    e = diag(D);
    r = sum(e >= options.theta * e(1));
    fprintf('Iter %d, obj:%0.8f, dinf:%0.1e, r:%d, p:%d, time:%0.2fs\n', ...
            iter, obj, dinf, r, p, toc(timespend));
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
        Y = V(:,1:r)' .* e(1:r);
        p = r;
    end
    tau = max(0.4*tau, 1e-8);
    nne = min(sum(dS < -tau), options.delta);
    if options.line_search == 1
        U = [zeros(p, n2); vS(:,1:nne)'];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y; zeros(nne, n2)];
    else
        Yn  = [Y; options.alpha * vS(:,1:nne)'];
        Y1n = Yn(:, 1:n);
        Y2n = Yn(:, n+1:n2);
        sc  = sqrt((sum(Y1n.^2) + sum(Y2n.^2)) / 2);
        Y   = [Y1n ./ sc, Y2n ./ sc];
    end
end

data.X        = X;
data.S        = S;
data.z        = z;
data.dinf     = dinf;
data.gradnorm = gradnorm;
data.time     = toc(timespend);
if data.status == 0 && dinf > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end
fprintf('ManiSDP_onlycomplexunitdiag: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    % --- cost function: f = 0.5*trace(C*X) = sum(eG_pair) ---
    function [f, store] = cost(Y, store)
        YC      = Y * C;
        Y1_     = Y(:, 1:n);
        Y2_     = Y(:, n+1:n2);
        eG1_    = sum(Y1_ .* YC(:, 1:n));
        eG2_    = sum(Y2_ .* YC(:, n+1:n2));
        eG_pair = (eG1_ + eG2_) / 2;   % updates shared closure
        f       = sum(eG_pair);
    end

    % --- Riemannian gradient (paired-oblique projection of Y*C) ---
    function [G, store] = grad(Y, store)
        G = [YC(:, 1:n)    - Y(:, 1:n)    .* eG_pair, ...
             YC(:, n+1:n2) - Y(:, n+1:n2) .* eG_pair];
    end

    % --- Riemannian Hessian ---
    % Hess f[U] = Proj(U*C) + Weingarten(U, rgrad f)
    % Proj_paired(W)_i  = W_i - [Y1_i;Y2_i] * (Y1_i'*W1_i + Y2_i'*W2_i)/2
    % Weingarten term   = -U .* eG_pair  (same structure as standard oblique)
    function [H, store] = hess(Y, U, store)
        eH     = U * C;
        eH1    = eH(:, 1:n);
        eH2    = eH(:, n+1:n2);
        Y1_    = Y(:, 1:n);
        Y2_    = Y(:, n+1:n2);
        lam_H  = (sum(Y1_.*eH1) + sum(Y2_.*eH2)) / 2;   % 1 x n
        H = [(eH1 - Y1_.*lam_H) - U(:, 1:n)    .* eG_pair, ...
             (eH2 - Y2_.*lam_H) - U(:, n+1:n2) .* eG_pair];
    end

    % --- line search along U from Y, with paired retraction ---
    function nY = line_search(Y, U)
        % The staircase direction satisfies Y'*U = 0. The derivative
        % in theta = alpha^2 is half the Hessian curvature of co along U;
        % co is twice the objective used by the Manopt callbacks.
        nY = Y;
        ls_pairnorm = sqrt((sum(Y(:,1:n).^2, 1) ...
            + sum(Y(:,n+1:n2).^2, 1))/2);
        if any(~isfinite(ls_pairnorm)) || any(ls_pairnorm == 0)
            return;
        end
        % Restore paired squared norms equal to two after rank truncation.
        Y = Y./[ls_pairnorm, ls_pairnorm];
        nY = Y;
        ls_pair_u2 = sum(U(:,1:n).^2, 1) ...
            + sum(U(:,n+1:n2).^2, 1);
        if any(~isfinite(ls_pair_u2)) || max(ls_pair_u2) == 0
            return;
        end
        ls_pairrate = ls_pair_u2/2;

        % Only equal diagonal coefficients within a pair are constant.
        % Removing them improves accuracy and preserves sparsity of C.
        ls_diagonal = diag(C);
        ls_pairdiag = 0.5*ls_diagonal(1:n) + 0.5*ls_diagonal(n+1:n2);
        ls_C = C - spdiags([ls_pairdiag; ls_pairdiag], 0, n2, n2);
        ls_YC = Y*ls_C;
        ls_columncost = sum(Y.*ls_YC, 1);
        ls_paircost = (ls_columncost(1:n) + ls_columncost(n+1:n2))/2;
        ls_slope = sum(U.*(U*ls_C), 'all') - sum(ls_paircost.*ls_pair_u2);
        if ~isfinite(ls_slope) || ls_slope >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(max(ls_pairrate)));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            ls_pairscale = sqrt(1 + ls_theta*ls_pairrate);
            ls_pairshift = (ls_theta*ls_pairrate)./(1 + ls_pairscale);
            ls_scale = [ls_pairscale, ls_pairscale];
            ls_shift = [ls_pairshift, ls_pairshift];
            % Stable paired-retraction displacement and exact cost change.
            ls_delta = (ls_alpha*U - Y.*ls_shift)./ls_scale;
            ls_deltaC = ls_delta*ls_C;
            ls_change = 2*sum(ls_delta.*ls_YC, 'all') ...
                + sum(ls_delta.*ls_deltaC, 'all');
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                ls_trial = Y + ls_alpha*U;
                ls_trialnorm = sqrt((sum(ls_trial(:,1:n).^2, 1) ...
                    + sum(ls_trial(:,n+1:n2).^2, 1))/2);
                nY = ls_trial./[ls_trialnorm, ls_trialnorm];
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep the normalized base point if no trial is acceptable.
    end

    function nY = retr_paired(Y, V)
        Y1_ = Y(:, 1:n)    + V(:, 1:n);
        Y2_ = Y(:, n+1:n2) + V(:, n+1:n2);
        sc  = sqrt((sum(Y1_.^2) + sum(Y2_.^2)) / 2);
        nY  = [Y1_ ./ sc, Y2_ ./ sc];
    end

    % --- paired-oblique manifold factory ---
    % Points: p x 2n matrices where ||col_i||^2 + ||col_{n+i}||^2 = 2 for i=1..n
    function M = pairedobliquefactory(p, n_)
        n2_ = 2 * n_;
        M.dim         = @() (2*p - 1) * n_;
        M.inner       = @(x, d1, d2) sum(d1 .* d2, 'all');
        M.norm        = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi * sqrt(n_);
        M.proj        = @proj_t;
        M.tangent     = @proj_t;
        M.retr        = @retr_m;
        M.rand        = @rand_pt;
        M.lincomb     = @matrixlincomb;
        M.zerovec     = @(x) zeros(p, n2_);
        M.transp      = @(x1, x2, d) proj_t(x2, d);

        function PU = proj_t(Y, U)
            lam = (sum(Y(:,1:n_)    .* U(:,1:n_)) + ...
                   sum(Y(:,n_+1:n2_).* U(:,n_+1:n2_))) / 2;   % 1 x n_
            PU  = [U(:,1:n_)     - Y(:,1:n_)     .* lam, ...
                   U(:,n_+1:n2_) - Y(:,n_+1:n2_) .* lam];
        end

        function y = retr_m(x, d)
            xtd1 = x(:,1:n_)     + d(:,1:n_);
            xtd2 = x(:,n_+1:n2_) + d(:,n_+1:n2_);
            sc   = sqrt((sum(xtd1.^2) + sum(xtd2.^2)) / 2);
            y    = [xtd1 ./ sc, xtd2 ./ sc];
        end

        function Yr = rand_pt()
            Y1_ = randn(p, n_);
            Y2_ = randn(p, n_);
            sc  = sqrt((sum(Y1_.^2) + sum(Y2_.^2)) / 2);
            Yr  = [Y1_ ./ sc, Y2_ ./ sc];
        end
    end
end
