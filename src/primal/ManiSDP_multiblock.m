% This function solves linear SDPs with unital diagonal.
% Min  <C, X>
% s.t. A(X) = b,
%      X in S_+^{n_1×...×n_t}
%      X_ii = 1, i = 1,...,n.

function [X, obj, data] = ManiSDP_multiblock(At, b, c, K, options)

n = K.s;
nb = length(n);
if ~isfield(options,'min_facsize'); options.min_facsize = 2; end
if ~isfield(options,'p0'); options.p0 = 2*ones(nb,1); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gama'); options.gama = 3; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e2; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-1; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tao'); options.tao = 10; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', max(n), size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');

A = At';
p = n;
for i = 1:nb
    if n(i) >= options.min_facsize
        p(i) = options.p0(i);
    end
end
sigma = options.sigma0;
gama = options.gama;
y = zeros(length(b), 1);
normb = 1 + norm(b);
x = zeros(sum(n.^2), 1);
YU = zeros(sum(n.^2), 1);
Y = [];
U = [];
X = cell(nb, 1);
S = cell(nb, 1);
% fac_size = [];
% seta = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgrad;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
%     fac_size = [fac_size; p];
    problem.M = multiblockmanifold(p, n, K.nob);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
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
%     if pinf >= gradnorm
        y = y - sigma*Axb;   
%     end
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
        S{i} = 0.5*(S{i}+S{i}');
        ind = ind + n(i)^2;
        [vS{i}, dS{i}] = eig(S{i}, 'vector');
        dinfs(i) = max(0, -dS{i}(1))/(1+abs(dS{i}(end)));
    end
    dinf = max(dinfs);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, p_max:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,   max(p),    sigma,       toc(timespend));
    eta = max([gap, pinf, dinf]);
%     seta = [seta; eta];
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 10) == 0
        if iter > 20 && gap > gap0 && pinf > pinf0 && dinf > dinf0
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
               nne = max(min(sum(dS{i} < 0), options.delta), 1);
            else
               nne = min(sum(dS{i} < 0), options.delta);
            end
            if p(i) + nne > n(i)
               nne = 0;
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
    if pinf < options.tao*gradnorm
        sigma = max(sigma/gama, options.sigma_min);
    else
        sigma = min(sigma*gama, options.sigma_max);
    end
%    tolgrad = pinf/2;
end
data.S = S;
data.y = y;
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

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

%     function Y = line_search(Y, U)
%         alpha = [0.0001;0.0002;0.0003;0.0004;0.0005];
%         val = zeros(length(alpha),1);
%         nY = cell(nb, 1);
%         for k = 1:length(alpha)
%             for i = 1:nb
%                 nY{i} = [Y{i}; alpha(k)*U{i}];
%                 if i <= K.nob
%                     nY{i} = nY{i}./sqrt(sum(nY{i}.^2));
%                 end
%             end
%             val(k) = co(nY);
%         end
%         [~, I] = min(val);
%         for i = 1:nb
%             Y{i} = [Y{i}; alpha(I)*U{i}];
%             if i <= K.nob
%                 Y{i} = Y{i}./sqrt(sum(Y{i}.^2));
%             end
%         end
%     end
        
    function val = co(Y)
        ind = 1;
        for i = 1:nb
            X{i} = Y{i}'*Y{i};
            x(ind:ind+n(i)^2-1) = X{i}(:);
            ind = ind + n(i)^2;
        end
        Axb = A*x - b - y/sigma;
        val = c'*x + 0.5*sigma*(Axb'*Axb);
    end


    function nY = line_search(Y, U)
         alpha = 0.5;
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
            store.eG{i} = sum(Y{i}.*G{i});
            if i <= K.nob
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
             H{i} = H{i} + 8*sigma*Y{i}*reshape(AyU(ind:ind+n(i)^2-1), n(i), n(i));
             if i <= K.nob
                 H{i} = H{i} - Y{i}.*sum(Y{i}.*H{i}) - U{i}.*store.eG{i};
             end
             ind = ind + n(i)^2;
        end
    end
end