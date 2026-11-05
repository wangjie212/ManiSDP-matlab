% This function solves SDPs with multi-blocks using the dual approach, where the first K.nob dual blocks (diag(S_i) = 1) are unit-diagonal:
%    sup  <C, X> + <c, w>
%    s.t. A(X) + B(w) = b
%         X in S_+^{n_1×...×n_t}
%         w in R^l
%         diag(S_i) = 1, i = 1,...,K.nob

function [X, obj, data] = ManiDSDP_multiblock(A, b, c, K, options)

if nargin < 5 || isempty(options)
    options = struct();
end

n = K.s;
nb = length(n);
if ~isfield(options,'min_facsize'); options.min_facsize = 2; end
if ~isfield(options,'p0'); options.p0 = ones(nb, 1); end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter = 1000; end
if ~isfield(options,'gamma'); options.gamma = 2; end
if ~isfield(options,'tau'); options.tau = 1e-2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.2; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'kappa1'); options.kappa1 = 1e1; end
if ~isfield(options,'kappa2'); options.kappa2 = 1e1; end
if ~isfield(options,'line_search'); options.line_search = 1; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', max(n), size(b,1));
warning_state = warning('off', 'manopt:trs_tCG_cached:memory');
warning_cleanup = onCleanup(@() warning(warning_state));

normc = 1 + norm(c);
if isfield(K, 'f')
    nf = K.f;
    B = A(:, 1:nf);
    cf = c(1:nf);
else
    nf = 0; 
end
A = A(:, nf+1:end);
c = c(nf+1:end);

if ~isfield(options, 'dAAt')
    iA = ((A*A')\A)'; % (inv(A*A')*A)'
else
    iA = (sparse(1:size(b,1),1:size(b,1),options.dAAt)\A)';
end
bA = iA*b;
if nf ~= 0; iAB = iA*B; end 

p = n;
for i = 1:nb
    if n(i) >= options.min_facsize
        p(i) = options.p0(i);
    end
end
sigma = options.sigma0; 
gamma = options.gamma;
tau = options.tau;
x = zeros(sum(n.^2), 1);
s = zeros(sum(n.^2), 1);
YU = zeros(sum(n.^2), 1);
if nf ~= 0; w = zeros(nf, 1); end
Y = [];
U = [];
S = cell(nb,1);
X = cell(nb,1);

problem.costgrad = @costgrad;
% problem.cost = @cost;
% problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = options.TR_maxinner;    % maximum Hessian calls per iteration
opts.maxiter   = options.TR_maxiter;
% opts.tolgradnorm = options.tolgradnorm;
% opts.trscache = false;
data.status = 0;

timespend = tic;
for iter = 1:options.ADMM_maxiter
    problem.M = multiblockmanifold(p, n, K.nob);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tau;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    ind = 1;
    for i = 1:nb
        S{i} = Y{i}'*Y{i};
        s(ind:ind+n(i)^2-1) = S{i}(:);
        ind = ind + n(i)^2;
    end
    sc = s - c;
    y = iA'*sc;
    As = A'*y - sc;
    nA = norm(As);
    if nf ~= 0
        Af = B'*y - cf;
        nA = nA + norm(Af);
    end
    pinf = nA/normc;
    by = b'*y;
    if K.nob == nb
        x = x - sigma*As;
    else
        x = x + sigma*(iAB*(Af-w/sigma) + A'*(iA'*(As-x/sigma)) - As);
    end
    if nf ~= 0
        w = w - sigma*Af;
        obj = c'*x + cf'*w;
    end
    dinfs = zeros(nb, 1);
    ind = 1;
    for i = 1:nb
        X{i} = reshape(x(ind:ind+n(i)^2-1)+bA(ind:ind+n(i)^2-1), n(i), n(i));
        if i <= K.nob
            z = sum(S{i}.*X{i});
            obj = obj + sum(z);
            X{i} = X{i} - diag(z);
        end
        ind = ind + n(i)^2;
        [vX{i}, dX{i}] = eig(X{i}, 'vector');
        dinfs(i) = max(0, -dX{i}(1))/(1+abs(dX{i}(end)));
    end
    dinf = max(dinfs);
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, p_max:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,    max(p),   sigma, toc(timespend));
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
            [~, D, V] = svd(Y{i}, 'econ');
            if size(D, 1) > 1
               e = diag(D);
            else
               e = D(1);
            end
            r = sum(e > options.theta*e(1));
            if r < p(i)
                Y{i} = diag(e(1:r))*V(:,1:r)';
                p(i) = r;
            end
            if p(i) == 1 && i <= K.nob
                nne = 1;
            else
                tau = max(0.4*tau, 1e-8);
                nne = min(sum(dX{i} < -tau), options.delta);
            end
            if options.line_search == 1
                U{i} = [zeros(p(i), n(i)); vX{i}(:,1:nne)'];
            end
            p(i) = p(i) + nne;
            if options.line_search == 1
               Y{i}  = [Y{i}; zeros(nne, n(i))];
            else
               Y{i} = [Y{i}; options.alpha*vX{i}(:,1:nne)'];
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

    function nY = line_search(Y, U)
        % Each staircase direction satisfies Y{i}'*U{i} = 0. The
        % pullback is even in alpha, and its derivative in theta =
        % alpha^2 is half the Hessian curvature on the product manifold.
        nY = Y;
        ls_base = Y;
        ls_U = cell(nb, 1);
        ls_S0 = cell(nb, 1);
        ls_UtU = cell(nb, 1);
        ls_u2 = cell(nb, 1);
        ls_indices = cell(nb, 1);
        ls_diagonal = cell(nb, 1);
        ls_s0 = zeros(sum(n.^2), 1);
        ls_ds0 = ls_s0;
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
            ls_S0{ls_block} = ls_base{ls_block}'*ls_base{ls_block};
            ls_UtU{ls_block} = ls_U{ls_block}'*ls_U{ls_block};
            ls_dS = ls_UtU{ls_block};
            if ls_block <= K.nob
                ls_diagonal{ls_block} = 1:(n(ls_block)+1):n(ls_block)^2;
                ls_S0{ls_block}(ls_diagonal{ls_block}) = 1;
                ls_dS = ls_dS - 0.5*(ls_u2{ls_block}.'.*ls_S0{ls_block} ...
                    + ls_S0{ls_block}.*ls_u2{ls_block});
                ls_dS(ls_diagonal{ls_block}) = 0;
                ls_step_scale = max(ls_step_scale, max(ls_u2{ls_block}));
            else
                ls_step_scale = max(ls_step_scale, sum(ls_u2{ls_block}));
            end
            ls_indices{ls_block} = ls_offset + (1:n(ls_block)^2);
            ls_s0(ls_indices{ls_block}) = ls_S0{ls_block}(:);
            ls_ds0(ls_indices{ls_block}) = ls_dS(:);
            ls_offset = ls_offset + n(ls_block)^2;
        end
        nY = ls_base;
        if ~isfinite(ls_step_scale) || ls_step_scale == 0
            return;
        end
        % Local residuals preserve the Manopt callbacks' shared state.
        ls_sc0 = ls_s0 - c;
        ls_y0 = iA'*ls_sc0;
        ls_As0 = A'*ls_y0 - ls_sc0 - x/sigma;
        ls_linear = iA*(b + sigma*(A*ls_As0)) - sigma*ls_As0;
        if nf ~= 0
            ls_Af0 = B'*ls_y0 - cf - w/sigma;
            ls_linear = ls_linear + sigma*(iA*(B*ls_Af0));
        end
        ls_slope = ls_linear'*ls_ds0;
        if ~isfinite(ls_slope) || ls_slope >= 0
            return;
        end

        ls_alpha = min(1, 1/sqrt(ls_step_scale));
        ls_armijo = 0.1;
        for ls_iter = 1:50
            ls_theta = ls_alpha^2;
            ls_ds = zeros(size(ls_s0));
            for ls_block = 1:nb
                ls_dS = ls_theta*ls_UtU{ls_block};
                if ls_block <= K.nob
                    ls_scale = sqrt(1 + ls_theta*ls_u2{ls_block});
                    ls_shift = (ls_theta*ls_u2{ls_block})./(1 + ls_scale);
                    % Exact normalized Gram change, without cancellation.
                    ls_dS = (ls_dS ...
                        - (ls_shift.' + ls_shift + ls_shift.'*ls_shift) ...
                          .*ls_S0{ls_block})./(ls_scale.'*ls_scale);
                    ls_dS(ls_diagonal{ls_block}) = 0;
                end
                ls_ds(ls_indices{ls_block}) = ls_dS(:);
            end
            % Aggregate before the quadratic term to include constraints
            % that couple different blocks.
            ls_dy = iA'*ls_ds;
            ls_dAs = A'*ls_dy - ls_ds;
            ls_quadratic = ls_dAs'*ls_dAs;
            if nf ~= 0
                ls_dAf = B'*ls_dy;
                ls_quadratic = ls_quadratic + ls_dAf'*ls_dAf;
            end
            ls_change = ls_linear'*ls_ds + 0.5*sigma*ls_quadratic;
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

    function [f, G, store] = costgrad(Y, store)
        ind = 1;
        for i = 1:nb
            S{i} = Y{i}'*Y{i};
            s(ind:ind+n(i)^2-1) = S{i}(:);
            ind = ind + n(i)^2;
        end
        sc = s - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        if nf ~= 0
           Af = B'*y - cf - w/sigma;
           f = b'*y + 0.5*sigma*(As'*As + Af'*Af);
        end
        if K.nob == nb
            tt = bA - sigma*As;
        else
            tt = bA + sigma*(iAB*Af + A'*(iA'*As) - As);
        end
        ind = 1;
        for i = 1:nb
            X{i} = reshape(tt(ind:ind+n(i)^2-1), n(i), n(i));
            G{i} = 2*Y{i}*X{i};
            if i <= K.nob
                store.eG{i} = G{i};
                G{i} = G{i} - Y{i}.*sum(Y{i}.*G{i});
            end
            ind = ind + n(i)^2;
        end
    end

    function [H, store] = hess(Y, U, store)
        ind = 1;
        for i = 1:nb
            T = U{i}'*Y{i};
            YU(ind:ind+n(i)^2-1) = T(:);
            ind = ind + n(i)^2;
            H{i} = 2*U{i}*X{i} + 2*sigma*Y{i}*(T + T');
        end
        if K.nob == nb
            tYU = -2*(YU'*iA)*A;
        else
            yAU = (A'*(iA'*YU))';
            tYU = -4*yAU + 2*(iAB*(iAB'*YU))' + 2*(A'*(iA'*yAU'))';
        end
        ind = 1;
        for i = 1:nb
             H{i} = H{i} + 2*sigma*Y{i}*reshape(tYU(ind:ind+n(i)^2-1), n(i), n(i));
             if i <= K.nob
                 H{i} = H{i} - Y{i}.*sum(Y{i}.*H{i}) - U{i}.*sum(Y{i}.*store.eG{i});
             end
             ind = ind + n(i)^2;
        end
    end

end
