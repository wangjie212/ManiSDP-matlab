% Solves complex unit-diagonal SDP via paired-oblique manifold on real conversion.
%
% Problem (real, after convertCtoR):
%   min  <C_r, X_r>
%   s.t. X_r >= 0
%        X_r(i,i) + X_r(n+i,n+i) = 2,  i = 1,...,n
%
% Factorization: X_r = Y_r'*Y_r, where Y_r is (p x 2n).
% Constraint encoded implicitly: ||col_i(Y_r)||^2 + ||col_{n+i}(Y_r)||^2 = 2.
%
% Analogous to ManiSDP_onlyunitdiag for real MaxCut; no explicit At matrix needed.
% Input C_r must be the 2n x 2n real matrix from convertCtoR (already divided by 2).

function [X, obj, data] = ManiSDP_onlycomplexunitdiag(C_r, options)

n2 = size(C_r, 1);   % = 2n
n  = n2 / 2;

if ~isfield(options,'p0');          options.p0          = 2;    end
if ~isfield(options,'AL_maxiter');  options.AL_maxiter  = 20;   end
if ~isfield(options,'tol');         options.tol         = 1e-8; end
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
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;

    % Compute objective and dual certificate
    X        = Y'*Y;
    Y1       = Y(:, 1:n);
    Y2       = Y(:, n+1:n2);
    YC_out   = Y * C_r;
    eG1      = sum(Y1 .* YC_out(:, 1:n));     % 1 x n
    eG2      = sum(Y2 .* YC_out(:, n+1:n2));  % 1 x n
    eG_pair  = (eG1 + eG2) / 2;               % 1 x n, dual multipliers y_i
    obj      = 2 * sum(eG_pair);              % = trace(C_r * X_r)

    % Dual slack: S_r = C_r - diag([y; y]), where y_i = 2*eG_pair_i
    % (from KKT: (Y*C_r)[:,i] = y_i/2 * Y[:,i], so y_i = 2*eG_pair_i,
    %  and S_r = C_r - sum_i y_i*(e_i*e_i'+e_{n+i}*e_{n+i}')/2
    %          = C_r - diag([eG_pair, eG_pair]))
    z  = [eG_pair, eG_pair];          % 1 x 2n
    S  = C_r - diag(z(:));            % 2n x 2n
    [vS, dS] = eig(full(S), 'vector');
    dS = sort(real(dS));
    dinf = max(0, -dS(1)) / (1 + abs(dS(end)));

    [~, D, V] = svd(Y);
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
    nne = max(min(sum(dS < 0), options.delta), 1);
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

    % --- cost function: f = 0.5*trace(C_r*X_r) = sum(eG_pair) ---
    function [f, store] = cost(Y, store)
        YC      = Y * C_r;
        Y1_     = Y(:, 1:n);
        Y2_     = Y(:, n+1:n2);
        eG1_    = sum(Y1_ .* YC(:, 1:n));
        eG2_    = sum(Y2_ .* YC(:, n+1:n2));
        eG_pair = (eG1_ + eG2_) / 2;   % updates shared closure
        f       = sum(eG_pair);
    end

    % --- Riemannian gradient (paired-oblique projection of Y*C_r) ---
    function [G, store] = grad(Y, store)
        G = [YC(:, 1:n)    - Y(:, 1:n)    .* eG_pair, ...
             YC(:, n+1:n2) - Y(:, n+1:n2) .* eG_pair];
    end

    % --- Riemannian Hessian ---
    % Hess f[U] = Proj(U*C_r) + Weingarten(U, rgrad f)
    % Proj_paired(W)_i  = W_i - [Y1_i;Y2_i] * (Y1_i'*W1_i + Y2_i'*W2_i)/2
    % Weingarten term   = -U .* eG_pair  (same structure as standard oblique)
    function [H, store] = hess(Y, U, store)
        eH     = U * C_r;
        eH1    = eH(:, 1:n);
        eH2    = eH(:, n+1:n2);
        Y1_    = Y(:, 1:n);
        Y2_    = Y(:, n+1:n2);
        lam_H  = (sum(Y1_.*eH1) + sum(Y2_.*eH2)) / 2;   % 1 x n
        H = [(eH1 - Y1_.*lam_H) - U(:, 1:n)    .* eG_pair, ...
             (eH2 - Y2_.*lam_H) - U(:, n+1:n2) .* eG_pair];
    end

    % --- co: unnormalized cost, used only in line search ---
    function val = co(Y)
        val = sum((Y * C_r) .* Y, 'all');
    end

    % --- line search along U from Y, with paired retraction ---
    function nY = line_search(Y, U)
        alpha = 1;
        cost0 = co(Y);
        i     = 1;
        nY    = retr_paired(Y, alpha * U);
        while i <= 15 && co(nY) - cost0 > -1e-3
            alpha = 0.8 * alpha;
            nY    = retr_paired(Y, alpha * U);
            i     = i + 1;
        end
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
