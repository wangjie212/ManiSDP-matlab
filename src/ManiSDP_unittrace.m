% This function solves linear SDPs with unital trace.
% Min  <C, X>
% s.t. A(X) = b,
%      tr(X) = 1.

function [Y, S, y, fval, error] = ManiSDP_unittrace(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 1;
sigma = 1e1;
sigma_min = 1e2; 
sigma_max = 1e7;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-8;
y = zeros(length(b),1);
normb = 1 + norm(b);
Y = [];
U = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;
opts.maxinner = 40;
opts.maxiter = 3;
timespend = tic;
for iter = 1:MaxIter
    problem.M = spherefactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tolgrad;
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y*Y';
    x = X(:);
    fval = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    eS = C - reshape(y'*A, n, n);
    lambda = sum(eS.*X,'all');
    S = eS - lambda*eye(n);
    [vS, dS] = eig(S, 'vector');
    mS = abs(min(dS))/(1+dS(end));
    by = b'*y + lambda;
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > 1e-2*e(1));
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       pinf,   gradnorm,    r,    p,    sigma,   toc(timespend));
    error = max([pinf, gap, mS]);
    if error < tao
        break;
    end
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = min(sum(dS < 0), 10);
    U = [zeros(n, p) vS(:,1:nne)];
    p = p + nne;
    Y = [Y zeros(n, nne)];
    % Y = [Y 0.04*vS(:,1:nne)];
    % Y = Y/norm(Y, 'fro');
%     if iter == 1 || neta > 0.5*eta
%         if sigma < 1e5
%               sigma = gama*sigma;
%         else
%               sigma = 1e3;
%         end
%     end
%     eta = neta;
    
    if pinf < gradnorm/6e3
          sigma = max(sigma/gama, sigma_min);
    else
          sigma = min(sigma*gama, sigma_max);
    end
%    tolgrad = 6e2*pinf;
end

    function val = co(Y)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
    end

    function nY = line_search(Y, U)
         alpha = 0.2;
         cost0 = co(Y);
         i = 1;
         nY = Y + alpha*U;
         nY = nY/norm(nY, 'fro');
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              nY = nY/norm(nY, 'fro');
              i = i + 1;
         end
    end

    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = c'*x + sigma/2*(Axb'*Axb);
        eS = C + sigma*reshape(Axb'*A, n, n);
        store.lambda = sum(X.*eS,'all');
        store.G = 2*eS*Y - 2*store.lambda*Y;
        store.eS = eS;
    end
    
    function [G, store] = grad(~, store)
        G = store.G;
    end

    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU= reshape(YU(:)'*At*A, n, n);
        H = 2*store.eS*U + 4*sigma*(AyU*Y);
        H = H - trace(H*Y')*Y - 2*store.lambda*U;
    end
end