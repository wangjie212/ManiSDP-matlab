function  [Y, fval, info] = SDP_ALM_subprog(At, b, C, c, n, p, sigma, y, Y0)
    % Pick the manifold of n-by-p matrices with unit norm rows.
    % manifold = obliquefactory(p, n, true);
    % manifold = elliptopefactory(n, p);
    manifold = symfixedrankYYfactory(n, p);
    % manifold = spectrahedronfactory(n, p);
    problem.M = manifold;
    
    % Define the cost function to be /minimized/.
    problem.cost = @cost;
    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        cx = x'*c;
        Axb = At'*x - b + y/sigma;
        f = cx + sigma*(Axb'*Axb);
        AxbA = Axb'*At';
        yA = reshape(AxbA, n, n);
        S = C + 2 * sigma*yA;
        store.G = 2*S*Y;
        store.S = S;
    end

    % Define the Riemannian gradient.
    problem.egrad = @egrad;
    function [G, store] = egrad(Y, store)
        G = store.G;
    end

    % If you want to, you can specify the Riemannian Hessian as well.
    problem.ehess = @ehess;
    function [H, store] = ehess(Y, Ydot, store)
        S = store.S;
        H = 2*S*Ydot;
        Xdot = Y*Ydot';
        xdot = Xdot(:);
        AxbdotA = xdot'*At*At';
        yAdot = reshape(AxbdotA, n, n);
        H = H + 8*sigma*(yAdot*Y);
    end

    % Call your favorite solver.
    opts = struct();
    opts.verbosity = 0;      % Set to 0 for no output, 2 for normal output
    opts.maxinner = 30;     % maximum Hessian calls per iteration
    opts.tolgradnorm = 1e-2; % tolerance on gradient norm
    opts.maxiter = 50;
    [Y, fval, info] = trustregions(problem, Y0, opts);

end
