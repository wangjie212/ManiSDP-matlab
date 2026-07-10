% Solves real SDP from a complex unit-diagonal SDP via convertCtoR:
%   min  <C_r, X_r>
%   s.t. A_r(X_r) = b_r
%        X_r >= 0
%        X_r(j,j) + X_r(n+j,n+j) = 2,  j = 1,...,n   (paired diagonal)
%
% Prepare inputs with: [At_r, b_r, c_r, K_r] = convertCtoR(At_c, b_c, c_c, K_c)
% where the original complex SDP has diag(X_c) = 1.

function [X, obj, data] = ManiSDP_complexunitdiag(At_r, b_r, c_r, K_r, options)

n2 = K_r.s;   % = 2n
n  = n2 / 2;

if ~isfield(options,'p0');          options.p0          = 2;     end
if ~isfield(options,'AL_maxiter');  options.AL_maxiter  = 300;   end
if ~isfield(options,'gama');        options.gama        = 2;     end
if ~isfield(options,'sigma0');      options.sigma0      = 1e-3;  end
if ~isfield(options,'sigma_min');   options.sigma_min   = 1e-2;  end
if ~isfield(options,'sigma_max');   options.sigma_max   = 1e7;   end
if ~isfield(options,'tol');         options.tol         = 1e-8;  end
if ~isfield(options,'theta');       options.theta       = 1e-3;  end
if ~isfield(options,'delta');       options.delta       = 8;     end
if ~isfield(options,'alpha');       options.alpha       = 0.1;   end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8;  end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20;    end
if ~isfield(options,'TR_maxiter');  options.TR_maxiter  = 4;     end
if ~isfield(options,'tau1');        options.tau1        = 1;     end
if ~isfield(options,'tau2');        options.tau2        = 1;     end
if ~isfield(options,'line_search'); options.line_search = 0;     end

fprintf('ManiSDP_complexunitdiag is starting...\n');
fprintf('Real SDP: n2 = %d (complex n = %d), m = %d\n', n2, n, length(b_r));
warning('off', 'manopt:trs_tCG_cached:memory');

A_r    = At_r';
p      = options.p0;
sigma  = options.sigma0;
gama   = options.gama;
y      = zeros(length(b_r), 1);
normb  = 1 + norm(b_r);
Y      = [];
U      = [];
Axb    = [];   % shared: cost -> grad
eS_    = [];   % shared: grad -> hess
fac_size = [];

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
    fac_size = [fac_size; p];
    problem.M = pairedobliquefactory(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;

    X    = Y'*Y;
    x    = X(:);
    obj  = c_r'*x;
    Axb  = A_r*x - b_r;
    pinf = norm(Axb) / normb;
    y    = y - sigma*Axb;

    % Dual certificate for paired-diagonal constraint
    eS_r = reshape(c_r - At_r*y, n2, n2);
    eS_r = (eS_r + eS_r') / 2;
    diag_XeS = sum(X .* eS_r, 1);                          % vectorized: O(n2²) vs O(p·n2²) for Y*eS_r
    z_pair = (diag_XeS(1:n) + diag_XeS(n+1:n2)) / 2;
    z_full = [z_pair, z_pair];
    S    = eS_r - diag(z_full(:));
    [vS, dS] = eig(full(S), 'vector');
    dS   = sort(real(dS));
    dinf = max(0, -dS(1)) / (1 + abs(dS(end)));

    by  = b_r'*y + sum(z_full);
    gap = abs(obj - by) / (abs(obj) + abs(by) + 1);
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e >= options.theta*e(1));

    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
            iter, obj, gap, pinf, dinf, gradnorm, r, p, sigma, toc(timespend));

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
            gap0 = gap; pinf0 = pinf; dinf0 = dinf;
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
    if pinf < options.tau1*gradnorm
        sigma = max(sigma/gama, options.sigma_min);
    elseif pinf > options.tau2*gradnorm
        sigma = min(sigma*gama, options.sigma_max);
    end
end

data.X        = X;
data.y        = y;
data.S        = S;
data.z        = z_full;
data.gap      = gap;
data.pinf     = pinf;
data.dinf     = dinf;
data.gradnorm = gradnorm;
data.time     = toc(timespend);
data.fac_size = fac_size;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end
fprintf('ManiSDP_complexunitdiag: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    % --- AL cost (Axb is shared with grad via nested-function workspace) ---
    function [f, store] = cost(Y, store)
        X_  = Y'*Y;
        x_  = X_(:);
        Axb = A_r*x_ - b_r - y/sigma;
        f   = c_r'*x_ + 0.5*sigma*(Axb'*Axb);
    end

    % --- Riemannian gradient on paired-oblique manifold ---
    function [G, store] = grad(Y, store)
        eS_ = reshape(c_r + sigma*At_r*Axb, n2, n2);
        eS_ = (eS_ + eS_') / 2;
        Y1_ = Y(:, 1:n);
        Y2_ = Y(:, n+1:n2);
        eG_ = 2*(Y * eS_);
        eG1 = eG_(:, 1:n);
        eG2 = eG_(:, n+1:n2);
        lam = (sum(Y1_.*eG1) + sum(Y2_.*eG2)) / 2;   % 1×n paired multipliers
        store.lam = lam;
        G = [(eG1 - Y1_.*lam), (eG2 - Y2_.*lam)];
    end

    % --- Riemannian Hessian on paired-oblique manifold ---
    function [H, store] = hess(Y, U, store)
        YU   = Y'*U;
        AyU  = reshape(At_r*(At_r'*YU(:)), n2, n2);
        eH_  = 2*U*eS_ + 4*sigma*(Y*AyU);          % uses shared eS_ from grad
        Y1_  = Y(:, 1:n);   Y2_  = Y(:, n+1:n2);
        U1_  = U(:, 1:n);   U2_  = U(:, n+1:n2);
        eH1  = eH_(:, 1:n); eH2  = eH_(:, n+1:n2);
        lamH = (sum(Y1_.*eH1) + sum(Y2_.*eH2)) / 2;
        lam  = store.lam;
        H = [(eH1 - Y1_.*lamH) - U1_.*lam, ...
             (eH2 - Y2_.*lamH) - U2_.*lam];
    end

    % --- AL cost for line search ---
    function val = co(Y)
        X_   = Y'*Y;
        x_   = X_(:);
        Axb_ = A_r*x_ - b_r - y/sigma;
        val  = c_r'*x_ + sigma/2*(Axb_'*Axb_);
    end

    % --- Line search along U with paired retraction ---
    function nY = line_search(Y, U)
        alpha = 1;
        cost0 = co(Y);
        i     = 1;
        nY    = retr_paired(Y, alpha*U);
        while i <= 15 && co(nY) - cost0 > -1e-3
            alpha = 0.8*alpha;
            nY    = retr_paired(Y, alpha*U);
            i     = i + 1;
        end
    end

    function nY = retr_paired(Y, V)
        Y1_ = Y(:,1:n)    + V(:,1:n);
        Y2_ = Y(:,n+1:n2) + V(:,n+1:n2);
        sc  = sqrt((sum(Y1_.^2) + sum(Y2_.^2)) / 2);
        nY  = [Y1_./sc, Y2_./sc];
    end

    % --- Paired-oblique manifold: ||col_j(Y)||^2 + ||col_{n+j}(Y)||^2 = 2 ---
    function M = pairedobliquefactory(p_, n_)
        n2_ = 2*n_;
        M.dim         = @() (2*p_ - 1)*n_;
        M.inner       = @(x, d1, d2) sum(d1 .* d2, 'all');
        M.norm        = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi*sqrt(n_);
        M.proj        = @proj_t;
        M.tangent     = @proj_t;
        M.retr        = @retr_m;
        M.rand        = @rand_pt;
        M.lincomb     = @matrixlincomb;
        M.zerovec     = @(x) zeros(p_, n2_);
        M.transp      = @(x1, x2, d) proj_t(x2, d);

        function PU = proj_t(Y_, U_)
            lam_ = (sum(Y_(:,1:n_)     .* U_(:,1:n_)) + ...
                    sum(Y_(:,n_+1:n2_) .* U_(:,n_+1:n2_))) / 2;
            PU = [U_(:,1:n_)     - Y_(:,1:n_)     .* lam_, ...
                  U_(:,n_+1:n2_) - Y_(:,n_+1:n2_) .* lam_];
        end

        function y_ = retr_m(x_, d_)
            xtd1 = x_(:,1:n_)     + d_(:,1:n_);
            xtd2 = x_(:,n_+1:n2_) + d_(:,n_+1:n2_);
            sc_  = sqrt((sum(xtd1.^2) + sum(xtd2.^2)) / 2);
            y_   = [xtd1./sc_, xtd2./sc_];
        end

        function Yr = rand_pt()
            Y1_ = randn(p_, n_);
            Y2_ = randn(p_, n_);
            sc_ = sqrt((sum(Y1_.^2) + sum(Y2_.^2)) / 2);
            Yr  = [Y1_./sc_, Y2_./sc_];
        end
    end
end
