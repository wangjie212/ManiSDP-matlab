% This function solves general linear SDPs.
% Min  <C, X>
% s.t. A(X) = b.

function [Y, S, y, fval, error] = ManiSDP(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 1;
sigma = 1e-2;
sigma_min = 1e-1; 
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
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.maxiter = 4;
opts.tolgradnorm = tolgrad;
timespend = tic;
for iter = 1:MaxIter
    problem.M = euclideanfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y*Y';
    x = X(:);
    fval = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    S = C - reshape(y'*A, n, n);
    [vS, dS] = eig(S, 'vector');
    mS = abs(min(dS))/(1+dS(end));
    by = b'*y;
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > 1e-2*e(1));
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       pinf,   gradnorm,    r,    p,    sigma,   toc(timespend));
    error = max([pinf, gap, mS]);
    if error < tao
        break;
    end
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = min(sum(dS < 0), 6);
    U = [zeros(n, p) vS(:,1:nne)];
    p = p + nne;
    Y = [Y zeros(n, nne)];
    % Y = [Y 0.02*vS(:,1:nne)];
        
%     if iter == 1 || pinf > 0.3*opinf
%         if sigma > 1e3
%             sigma = 1e2;
%         else
%             sigma = sigma*gama;
%         end
%     end
%     opinf = pinf;
    
    if pinf < gradnorm/1e2
          sigma = max(sigma/gama, sigma_min);
    else
          sigma = min(sigma*gama, sigma_max);
    end
%    tolgrad = pinf;
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
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              i = i + 1;
         end
    end

   function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = c'*x + sigma/2*(Axb'*Axb);
    end
    
    function [G, store] = grad(Y, store)
        S = C + sigma*reshape(Axb'*A, n, n);
        G = 2*S*Y;
    end

    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU = reshape(YU(:)'*At*A, n, n);
        H = 2*S*U + 4*sigma*(AyU*Y);
    end
end