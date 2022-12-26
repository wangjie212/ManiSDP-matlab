function  [Y, dinf] = SDP_subprog_sphere(At, A, b, c, C, n, p, sigma, y, Y0, U, tolgrad)
    manifold = spherefactory(n, p);
    problem.M = manifold;
    
    problem.cost = @cost;
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
    
    problem.grad = @grad;
    function [G, store] = grad(~, store)
        G = store.G;
    end

    problem.hess = @hess;
    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU= reshape(YU(:)'*At*A, n, n);
        H = 2*store.eS*U + 4*sigma*(AyU*Y);
        H = H - trace(H*Y')*Y - 2*store.lambda*U;
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

    if ~isempty(U)
        Y0 = line_search(Y0, U);
    end
    opts.verbosity = 0;
    opts.maxinner = 40;
    opts.tolgradnorm = tolgrad;
    opts.maxiter = 3;
    [Y, ~, info] = trustregions(problem, Y0, opts);
    dinf = info(end).gradnorm;
end
