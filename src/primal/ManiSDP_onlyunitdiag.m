% This function solves linear SDPs with unital diagonal and without extra affine constraints.
% Min  <C, X>
% s.t. X >= 0,
%      X_ii = 1, i = 1,...,n.

function [X, obj, data] = ManiSDP_onlyunitdiag(C, options)

if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 20; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-1; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.5; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 100; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 40; end
if ~isfield(options,'line_search'); options.line_search = 0; end

fprintf('ManiSDP is starting...\n');
n = size(C,1);
fprintf('SDP size: n = %i, m = %i\n', n, n);

p = options.p0;
Y = [];
U = [];
YC = [];
eG = [];
% fac_size = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgrad; % tolerance on gradient norm

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
%     fac_size = [fac_size; p];
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y'*Y;
    z = sum(C.*X);
    % z = sum((Y*C).*Y);
    obj = full(sum(z));
    S = C - diag(z);
    [vS, dS] = eig(full(S), 'vector');
    dinf = max(0, -dS(1))/(1+dS(end));
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, dinf:%0.1e, r:%d, p:%d, time:%0.2fs\n', ...
             iter,    obj,       dinf,       r,    p,    toc(timespend));
    if dinf < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 10) == 0
        if iter > 20 && dinf > dinf0
            data.status = 2;
            fprintf('Slow progress!\n');
            break;
        else
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
end
data.S = S;
data.z = z;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
% data.fac_size = fac_size;
if data.status == 0 && dinf > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    function val = co(Y)
        val = sum((Y*C).*Y,'all');
    end
    
%    function Y = line_search(Y, U)
%         alpha = [0.02;0.04;0.06;0.08;0.1;0.2];
%         val = zeros(length(alpha),1);
%         for i = 1:length(alpha)
%             nY = Y + alpha(i)*U;
%             nY = nY./sqrt(sum(nY.^2));
%             val(i) = co(nY);
%         end
%         [~,I] = min(val);
%         Y = Y + alpha(I)*U;
%         Y = Y./sqrt(sum(Y.^2));
%    end

    function nY = line_search(Y, U)
         alpha = 0.5;
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
        YC = Y*C;
        eG = sum(YC.*Y);
        f = 0.5*sum(eG);
    end

    function [G, store] = grad(Y, store)
        G = YC - Y.*eG;
    end

    function [He, store] = hess(Y, U, store)
        H = U*C;
        He = H - Y.*sum(Y.*H) - U.*eG;
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) sum(d1.*d2,'all');
%        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d,'fro');
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