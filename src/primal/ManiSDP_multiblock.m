% This function solves SDPs with multi-blocks using the primal approach, where the first K.nob blocks are unit-diagonal:
% Min  <C, X>
% s.t. A(X) = b
%      X in S_+^{n_1×...×n_t}
%      diag(X_i) = 1, i = 1,...,K.nob

function [X, obj, data] = ManiSDP_multiblock(At, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end

n = K.s;
nb = length(n);
if ~isfield(options,'min_facsize'); options.min_facsize = 2; end
if ~isfield(options,'p0'); options.p0 = ones(nb,1); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 1000; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 1e1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1e1; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', max(n), size(b,1));
warning_state = warning('off', 'manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));

A = At';
p = n;
for i = 1:nb
    if n(i) >= options.min_facsize
        p(i) = options.p0(i);
    end
end
sigma = options.sigma0;
gamma = options.gamma;
tau = options.tau;
y = zeros(length(b), 1);
normb = 1 + norm(b);
x = zeros(sum(n.^2), 1);
YU = zeros(sum(n.^2), 1);
Y = [];
U = [];
X = cell(nb, 1);
S = cell(nb, 1);
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
% opts.tolgradnorm = options.tolgradnorm;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    problem.M = multiblockmanifold(p, n, K.nob);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    ind = 1;
    for i = 1:nb
        X{i} = Y{i}'*Y{i};
        x(ind:ind+n(i)^2-1) = X{i}(:);
        ind = ind + n(i)^2;
    end
    obj = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;   
    cy = c - At*y;
    by = b'*y;
    dinfs = zeros(nb, 1);
    ind = 1;
    for i = 1:nb
        S{i} = reshape(cy(ind:ind+n(i)^2-1), n(i), n(i));
        if i <= K.nob
            z = sum(X{i}.*S{i});
            by = by + sum(z);
            S{i} = S{i} - diag(z);
        end
        ind = ind + n(i)^2;
        [vS{i}, dS{i}] = eig(S{i}, 'vector');
        dinfs(i) = max(0, -dS{i}(1))/(1+abs(dS{i}(end)));
    end
    dinf = max(dinfs);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, p_max:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,   max(p),    sigma,       toc(timespend));
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
    for i = 1:nb
        if n(i) >= options.min_facsize
            if p(i) > 1
                [~, D, V] = svd(Y{i}, 'econ');
                if size(D, 1) > 1
                   e = diag(D);
                else
                   e = D(1);
                end
                r = sum(e >= options.theta*e(1));
                if r < p(i)
                   Y{i} = diag(e(1:r))*V(:,1:r)';
                   p(i) = r;
                end
            end
            if p(i) == 1 && i <= K.nob
                nne = 1;
            else
                tau = max(0.4*tau, 1e-8);
                nne = min(sum(dS{i} < -tau), options.delta);
            end
            if options.line_search == 1
               U{i} = [zeros(p(i), n(i)); vS{i}(:,1:nne)'];
            end
            p(i) = p(i) + nne;
            if options.line_search == 1
                Y{i} = [Y{i}; zeros(nne,n(i))];
            else
                Y{i} = [Y{i}; options.alpha*vS{i}(:,1:nne)'];
                if i <= K.nob
                   Y{i} = Y{i}./sqrt(sum(Y{i}.^2));
                end
            end
        end
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
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    function nY = line_search(Y, U)
        % Each staircase direction satisfies Y{i}'*U{i} = 0. The
        % pullback is even in alpha, and its derivative in theta =
        % alpha^2 is half the Hessian curvature on the product manifold.
        nY = Y;
        ls_base = Y;
        ls_U = cell(nb, 1);
        ls_X0 = cell(nb, 1);
        ls_UtU = cell(nb, 1);
        ls_u2 = cell(nb, 1);
        ls_indices = cell(nb, 1);
        ls_diagonal = cell(nb, 1);
        ls_x0 = zeros(sum(n.^2), 1);
        ls_dx0 = ls_x0;
        ls_offset = 0;
        ls_step_scale = 0;
        for ls_block = 1:nb
            if any(~isfinite(ls_base{ls_block}), 'all')
                return;
            end
            if ls_block <= K.nob
                ls_norm = sqrt(sum(ls_base{ls_block}.^2, 1));
                if any(~isfinite(ls_norm)) || any(ls_norm == 0)
                    return;
                end
                ls_base{ls_block} = ls_base{ls_block}./ls_norm;
            end
            % Blocks below min_facsize can have empty or missing U cells.
            if ls_block <= numel(U) && ~isempty(U{ls_block})
                ls_U{ls_block} = U{ls_block};
            else
                ls_U{ls_block} = zeros(size(ls_base{ls_block}));
            end
            ls_u2{ls_block} = sum(ls_U{ls_block}.^2, 1);
            if any(~isfinite(ls_u2{ls_block}))
                return;
            end
            ls_X0{ls_block} = ls_base{ls_block}'*ls_base{ls_block};
            ls_UtU{ls_block} = ls_U{ls_block}'*ls_U{ls_block};
            ls_dX = ls_UtU{ls_block};
            if ls_block <= K.nob
                ls_diagonal{ls_block} = 1:(n(ls_block)+1):n(ls_block)^2;
                ls_X0{ls_block}(ls_diagonal{ls_block}) = 1;
                ls_dX = ls_dX - 0.5*(ls_u2{ls_block}.'.*ls_X0{ls_block} ...
                    + ls_X0{ls_block}.*ls_u2{ls_block});
                ls_dX(ls_diagonal{ls_block}) = 0;
                ls_step_scale = max(ls_step_scale, max(ls_u2{ls_block}));
            else
                ls_step_scale = max(ls_step_scale, sum(ls_u2{ls_block}));
            end
            ls_indices{ls_block} = ls_offset + (1:n(ls_block)^2);
            ls_x0(ls_indices{ls_block}) = ls_X0{ls_block}(:);
            ls_dx0(ls_indices{ls_block}) = ls_dX(:);
            ls_offset = ls_offset + n(ls_block)^2;
        end
        nY = ls_base;
        if ~isfinite(ls_step_scale) || ls_step_scale == 0
            return;
        end
        ls_residual = A*ls_x0 - b - y/sigma;
        ls_linear = c + sigma*(A'*ls_residual);
        ls_slope = ls_linear'*ls_dx0;
        if ~isfinite(ls_slope) || ls_slope >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(ls_step_scale));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            ls_dx = zeros(size(ls_x0));
            for ls_block = 1:nb
                ls_dX = ls_theta*ls_UtU{ls_block};
                if ls_block <= K.nob
                    ls_scale = sqrt(1 + ls_theta*ls_u2{ls_block});
                    ls_shift = (ls_theta*ls_u2{ls_block})./(1 + ls_scale);
                    % Exact normalized Gram change, without cancellation.
                    ls_dX = (ls_dX ...
                        - (ls_shift.' + ls_shift + ls_shift.'*ls_shift) ...
                          .*ls_X0{ls_block})./(ls_scale.'*ls_scale);
                    ls_dX(ls_diagonal{ls_block}) = 0;
                end
                ls_dx(ls_indices{ls_block}) = ls_dX(:);
            end
            % Aggregate before the quadratic term to include constraints
            % that couple different blocks.
            ls_Adx = A*ls_dx;
            ls_change = ls_linear'*ls_dx + 0.5*sigma*(ls_Adx'*ls_Adx);
            if isfinite(ls_change) && ls_change < 0 ...
                    && ls_change <= ls_armijo*ls_theta*ls_slope
                for ls_block = 1:nb
                    % Y is already padded to the direction's row count.
                    nY{ls_block} = ls_base{ls_block} + ls_alpha*ls_U{ls_block};
                    if ls_block <= K.nob
                        nY{ls_block} = nY{ls_block} ...
                            ./sqrt(sum(nY{ls_block}.^2, 1));
                    end
                end
                return;
            end
            ls_alpha = 0.5*ls_alpha;
        end
        % Keep the feasible base point if no trial is acceptable.
    end

    function [f, store] = cost(Y, store)
        ind = 1;
        for i = 1:nb
            X{i} = Y{i}'*Y{i};
            x(ind:ind+n(i)^2-1) = X{i}(:);
            ind = ind + n(i)^2;
        end
        Axb = A*x - b - y/sigma;
        f = c'*x + 0.5*sigma*(Axb'*Axb);
    end

    function [G, store] = grad(Y, store)
        tt = c + sigma*At*Axb;
        ind = 1;
        for i = 1:nb
            S{i} = reshape(tt(ind:ind+n(i)^2-1), n(i), n(i));
            G{i} = 2*Y{i}*S{i};
            if i <= K.nob
                store.eG{i} = sum(Y{i}.*G{i});
                G{i} = G{i} - Y{i}.*store.eG{i};
            end
            ind = ind + n(i)^2;
        end
    end

    function [H, store] = hess(Y, U, store)
        ind = 1;
        for i = 1:nb
            T = Y{i}'*U{i};
            YU(ind:ind+n(i)^2-1) = T(:);
            ind = ind + n(i)^2;
            H{i} = 2*U{i}*S{i};
        end
        AyU = (YU'*At)*A;
        ind = 1;
        for i = 1:nb
             H{i} = H{i} + 4*sigma*Y{i}*reshape(AyU(ind:ind+n(i)^2-1), n(i), n(i));
             if i <= K.nob
                 H{i} = H{i} - Y{i}.*sum(Y{i}.*H{i}) - U{i}.*store.eG{i};
             end
             ind = ind + n(i)^2;
        end
    end
end