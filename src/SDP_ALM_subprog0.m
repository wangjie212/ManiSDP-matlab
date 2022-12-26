function  [Y, dinf] = SDP_ALM_subprog0(At, A, b, c, C, n, p, sigma, y, Y0, U, tolgrad)
    manifold = euclideanfactory(n, p);
    problem.M = manifold;
    
    problem.cost = @cost;
    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = c'*x + sigma/2*(Axb'*Axb);
        S = C + sigma*reshape(Axb'*A, n, n);
        store.G = 2*S*Y;
        store.S = S;
    end
    
    problem.grad = @grad;
    function [G, store] = grad(~, store)
        G = store.G;
    end

    problem.hess = @hess;
    function [H, store] = hess(Y, U, store)
        Xdot = U*Y';
        yAdot = reshape(Xdot(:)'*At*A, n, n);
        H = 2*store.S*U + 4*sigma*(yAdot*Y);
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
   
    if ~isempty(U)
        Y0 = line_search(Y0, U);
    end
    opts.verbosity = 0;
    opts.maxinner = 40;
    opts.tolgradnorm = tolgrad;
    opts.maxiter = 20;
    [Y, ~, info] = trustregions(problem, Y0, opts);
    dinf = info(end).gradnorm;
end
