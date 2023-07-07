% This function solves linear SDPs with unital diagonal.
% Min  <C, X>
% s.t. A(X) = b,
%      X in S_+^{n_1×...×n_t}
%      X_ii = 1, i = 1,...,n.

function [X, obj, data] = ManiSDP_unitdiag_multiblock(At, b, c, K, options)

n = K.s;
nb = length(n);
if ~isfield(options,'p0'); options.p0 = ones(nb,1); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-3; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-3; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tao'); options.tao = 1; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', max(n), size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');

A = At';
p = options.p0;
sigma = options.sigma0;
gama = options.gama;
y = zeros(length(b), 1);
normb = 1 + norm(b);
x = zeros(sum(n.^2), 1);
YU = zeros(sum(n.^2), 1);
Y = [];
U = [];
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
for i = 1:nb
    M.(['M' num2str(i)]) = obliquefactoryNTrans(p(i), n(i));
end
elems = fieldnames(M);
timespend = tic;
for iter = 1:options.AL_maxiter
%     fac_size = [fac_size; p];
    problem.M = productmanifold(M);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    ind = 1;
    for i = 1:nb
        X{i} = Y.(elems{i})'*Y.(elems{i});
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
        eS{i} = reshape(cy(ind:ind+n(i)^2-1), n(i), n(i));
        z = sum(X{i}.*eS{i});
        by = by + sum(z);
        S{i} = eS{i} - diag(z);
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
        [~, D, V] = svd(Y.(elems{i}));
        e = diag(D);
        r = sum(e > options.theta*e(1));
        if r <= p(i) - 1         
           Y.(elems{i}) = V(:,1:r)'.*e(1:r);
           p(i) = r;
        end
        nne = max(min(sum(dS{i} < 0), options.delta), 1);
        if options.line_search == 1
            U.(elems{i}) = [zeros(p(i), n(i)); vS{i}(:,1:nne)'];
        end
        p(i) = p(i) + nne;
        if options.line_search == 1
            Y.(elems{i}) = [Y.(elems{i}); zeros(nne,n(i))];
        else
            Y.(elems{i}) = [Y.(elems{i}); options.alpha*vS{i}(:,1:nne)'];
            Y.(elems{i}) = Y.(elems{i})./sqrt(sum(Y.(elems{i}).^2));
        end
        if pinf < options.tao*gradnorm
            sigma = max(sigma/gama, options.sigma_min);
        else
            sigma = min(sigma*gama, options.sigma_max);
        end
    end
%    tolgrad = pinf/2;
end
data.S = S;
data.y = y;
data.z = z;
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
%         alpha = [0.02;0.04;0.06;0.08;0.1];
%         val = zeros(length(alpha),1);
%         for i = 1:length(alpha)
%             nY = Y + alpha(i)*U;
%             nY = nY./sqrt(sum(nY.^2));
%             val(i) = co(nY);
%         end
%         [~,I] = min(val);
%         Y = Y + alpha(I)*U;
%         Y = Y./sqrt(sum(Y.^2));
%     end
        
    function val = co(Y)
        X = Y'*Y;
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
    end

    function nY = line_search(Y, U)
         alpha = 0.1;
         cost0 = co(Y);
         i = 1;
         nY = Y + alpha*U;
         nY = nY./sqrt(sum(nY.^2));
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              nY = nY./sqrt(sum(nY.^2));
              i = i + 1;
         end
    end

    function [f, store] = cost(Y, store)
        ind = 1;
        for i = 1:nb
            X{i} = Y.(elems{i})'*Y.(elems{i});
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
            eS{i} = reshape(tt(ind:ind+n(i)^2-1), n(i), n(i));
            eG = 2*Y.(elems{i})*eS{i};
            store.YeG{i} = sum(Y.(elems{i}).*eG);
            G.(elems{i}) = eG - Y.(elems{i}).*store.YeG{i};
            ind = ind + n(i)^2;
        end
    end

    function [H, store] = hess(Y, U, store)
        ind = 1;
        for i = 1:nb
            T = Y.(elems{i})'*U.(elems{i});
            YU(ind:ind+n(i)^2-1) = T(:);
            ind = ind + n(i)^2;
            eH{i} = 2*U.(elems{i})*eS{i};
        end
        AyU = YU'*At*A;
        ind = 1;
        for i = 1:nb
             eH{i} = eH{i} + 4*sigma*Y.(elems{i})*reshape(AyU(ind:ind+n(i)^2-1), n(i), n(i));
             H.(elems{i}) = eH{i} - Y.(elems{i}).*sum(Y.(elems{i}).*eH{i}) - U.(elems{i}).*store.YeG{i};
             ind = ind + n(i)^2;
        end
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
%        M.inner = @(x, d1, d2) sum(d1.*d2,'all');
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d, 'fro');
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U);

        M.retr = @retraction;
        function y = retraction(x, d, t)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2));
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d);

        function x = random(n, m)
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end