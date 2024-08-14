% This function solves SDPs with multi-blocks using the dual approach, where the first K.nob blocks are unit-diagonal:
%    sup  <C, X> + <c, x>
%    s.t. A(X) + B(x) = b
%         X in S_+^{n_1×...×n_t}
%         diag(X_i) = 1, i = 1,...,K.nob
%         x in R^l

function [X, obj, data] = ManiDSDP_multiblock(A, b, c, K, options)

n = K.s;
nb = length(n);
if ~isfield(options,'min_facsize'); options.min_facsize = 2; end
if ~isfield(options,'p0'); options.p0 = ones(nb,1); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 1000; end
if ~isfield(options,'gama'); options.gama = 3; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.2; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tao'); options.tao = 10; end
if ~isfield(options,'line_search'); options.line_search = 1; end

fprintf('ManiDSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', max(n), size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');
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

if ~isfield(options,'dAAt'); options.dAAt = diag(A*A'); end
iA = (sparse(1:size(b,1),1:size(b,1),options.dAAt)\A)';
% iA = ((A*A')\A)';
bA = iA*b;
if nf ~= 0; iAB = iA*B; end 

p = n;
for i = 1:nb
    if n(i) >= options.min_facsize
        p(i) = options.p0(i);
    end
end
sigma = options.sigma0; 
gama = options.gama;
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
opts.tolgradnorm = options.tolgrad;
% opts.trscache = false;
data.status = 0;

timespend = tic;
for iter = 1:options.AL_maxiter
    problem.M = multiblockmanifold(p, n, K.nob);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
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
        X{i} = 0.5*(X{i}+X{i}');
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
    for i = 1:nb
        if n(i) >= options.min_facsize
            [~, D, V] = svd(Y{i});
            e = diag(D);
            r = sum(e > options.theta*e(1));
            if r <= p(i) - 1
                Y{i} = V(:,1:r)'.*e(1:r);
                p(i) = r;
            end
            if i <= K.nob
                nne = max(min(sum(dX{i} < 0), options.delta), 1);
            else
                nne = min(sum(dX{i} < 0), options.delta);
            end
            if p(i) + nne > n(i)
               nne = 0;
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
    if pinf < options.tao*gradnorm
        sigma = max(sigma/gama, options.sigma_min);
    else
        sigma = min(sigma*gama, options.sigma_max);
    end
end
data.S = S;
data.y = y;
data.w = w;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
% data.fac_size = fac_size;
% data.seta = seta;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

        function val = co(Y)
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
               val = b'*y + 0.5*sigma*(As'*As + Af'*Af);
            end
        end


    function nY = line_search(Y, U)
         alpha = 1;
         cost0 = co(Y);
         nY = cell(nb, 1);
         for i = 1:nb
              nY{i} = [nY{i}; alpha*U{i}];
              if i <= K.nob
                  nY{i} = nY{i}./sqrt(sum(nY{i}.^2));
              end
         end
         k = 1;
         while k <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              for i = 1:nb
                  nY{i} = [Y{i}; alpha*U{i}];
                  if i <= K.nob
                      nY{i} = nY{i}./sqrt(sum(nY{i}.^2));
                  end
              end 
              k = k + 1;
         end
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
