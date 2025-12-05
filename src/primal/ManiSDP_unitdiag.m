% This function solves SDPs with unit-diagonal using the primal approach:
% Min  <C, X>
% s.t. A(X) = b
%      X >= 0
%      diag(X) = 1

function [X, obj, data] = ManiSDP_unitdiag(At, b, c, K, options)

n = K.s;
if ~isfield(options,'p0'); options.p0 = 2; end
% if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-3; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tau1'); options.tau1 = 1; end
if ~isfield(options,'tau2'); options.tau2 = 1; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');

A = At';
p = options.p0;
sigma = options.sigma0;
gama = options.gama;
y = zeros(length(b), 1);
normb = 1 + norm(b);
Y = [];
U = [];
fac_size = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
    fac_size = [fac_size; p];
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    x = X(:);
    obj = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;     
    y = y - sigma*Axb;
    eS = reshape(c - At*y, n, n);
    z = sum(X.*eS);
    S = eS - diag(z);
    [vS, dS] = eig(S, 'vector');
    dinf = max(0, -dS(1))/(1+dS(end));
    by = b'*y + sum(z);
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e >= options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,       toc(timespend));
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
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    nne = max(min(sum(dS < 0), options.delta), 1);
    if options.line_search == 1
       U = [zeros(p, n); vS(:,1:nne)'];
    end
    p = p + nne;
    if options.line_search == 1
       Y = [Y; zeros(nne,n)];
    else
       Y = [Y; options.alpha*vS(:,1:nne)'];
       Y = Y./sqrt(sum(Y.^2));
    end
    if pinf < options.tau1*gradnorm
          sigma = max(sigma/gama, options.sigma_min);
    elseif pinf > options.tau2*gradnorm
          sigma = min(sigma*gama, options.sigma_max);
    end
end
data.X = X;
data.y = y;
data.S = S;
data.z = z;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
data.fac_size = fac_size;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));  

    function val = co(Y)
        X = Y'*Y;
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
    end

    function nY = line_search(Y, U)
         alpha = 1;
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
        X = Y'*Y;
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = c'*x + 0.5*sigma*(Axb'*Axb);
    end

    function [G, store] = grad(Y, store)
        eS = reshape(c+sigma*At*Axb, n, n);
        eG = 2*Y*eS;
        store.YeG = sum(Y.*eG);
        G = eG - Y.*store.YeG;
    end

    function [H, store] = hess(Y, U, store)
        YU = Y'*U;
        AyU = reshape(A'*(At'*YU(:)), n, n);
        eH = 2*U*eS + 4*sigma*(Y*AyU);
        H = eH - Y.*sum(Y.*eH) - U.*store.YeG;
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
%        M.inner = @(x, d1, d2) sum(d1.*d2,'all');
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d, 'fro');
%        M.norm = @(x, d) norm(d);
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
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end