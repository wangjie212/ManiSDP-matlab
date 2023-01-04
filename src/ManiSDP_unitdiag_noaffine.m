% This function solves linear SDPs with unital diagonal and without extra affine constraints.
% Min  <C, X>
% s.t. X >= 0,
%      X_ii = 1, i = 1,...,n.

function [Y, S, fval, mS] = ManiSDP_unitdiag_noaffine(C)
n = size(C,1);
c = C(:);
p = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-8;
egrad = zeros(p, n);
Y = [];
U = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = 40;     % maximum Hessian calls per iteration
opts.tolgradnorm = tolgrad; % tolerance on gradient norm
opts.maxiter = 25;

timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    Y = trustregions(problem, Y, opts);
    X = Y'*Y;
    lambda = sum(C.*X);
    fval = full(sum(lambda));
    S = C - diag(lambda);
    [vS, dS] = eig(full(S), 'vector');
    mS = abs(min(dS))/(1+dS(end));
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > 1e-1*e(1));
    fprintf('Iter:%d, fval:%0.8f, mineigS:%0.1e, r:%d, p:%d, time:%0.2fs\n', ...
             iter,    fval,       mS,       r,    p,    toc(timespend));
    if mS < tao
        break;
    end
    if r <= p - 1         
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    nne = max(min(sum(dS < 0), 8), 1);
%    U = [zeros(p, n); vS(:,1:nne)'];
    p = p + nne;
%   Y = [Y; zeros(nne,n)]; 
    Y = [Y; 0.5*vS(:,1:nne)'];
    Y = Y./sqrt(sum(Y.^2));
end

    function val = co(Y)
        X = Y'*Y;
        val = c'*X(:);
    end

    function nY = line_search(Y, U)
         alpha = 0.2;
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
        f = c'*X(:);
    end

    function [G, store] = grad(Y, store)
        egrad = 2*Y*C;
        G = egrad - Y.*sum(Y.*egrad);
    end

    % If you want to, you can specify the Riemannian Hessian as well.
    function [He, store] = hess(Y, U, store)
        H = 2*U*C;
        He = H - Y.*sum(Y.*H) - U.*sum(Y.*egrad);
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d(:));
        M.tangent = @(X, U) U - X.*sum(X.*U);

        M.retr = @retraction;
        % Retraction on the oblique manifold
        function y = retraction(x, d)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2)); % y = normalize_columns(x + td);
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);

        % Uniform random sampling on the sphere.
        function x = random(n, m)
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end