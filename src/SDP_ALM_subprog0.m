function  Y = SDP_ALM_subprog0(At, A, b, c, C, n, p, sigma, y, Y0, tolgrad)
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
        S = store.S;
        H = 2*S*U;
        Xdot = U*Y';
        yAdot = reshape(Xdot(:)'*At*A, n, n);
        H = H + 4*sigma*(yAdot*Y);
    end

    opts.verbosity = 0;
    opts.maxinner = 20;
    opts.tolgradnorm = tolgrad;
    opts.maxiter = 5;
    Y = trustregions(problem, Y0, opts);
end
