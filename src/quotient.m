function  Y = quotient(At, A, b, c, C, n, p, sigma, y, Y0, U, tolgrad)  
    function M = quotientfactory(n, p)
        M.dim = @() n*p-0.5*p*(p-1);
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d(:));
        M.proj = @project;
        M.tangent = @project;
        M.retr = @retraction;
        M.rand = @() randn(n, p);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, p);
        
        function y = retraction(x, d)            
            y = x + d;
        end
        
        function P = project(Y, Z)
            X = Y'*Y;
            Omega = sylvester(X, X, Y'*Z-Z'*Y);
            P = Z - Y*Omega;
        end
    end

    problem.cost = @cost;
    function [f, store] = cost(Y, store)
        X0 = Y*Y';
        x0 = X0(:);
        Axb = A*x0 - b - y/sigma;
        store.Axb = Axb;
        f = c'*x0 + 0.5*sigma*(Axb'*Axb);
    end
    
    problem.grad = @grad;
    function [G, store] = grad(Y, store)
        yA0 = reshape(At*store.Axb, n, n);
        S = C + sigma*yA0;
        eG = 2*S*Y;
        G = problem.M.proj(Y, eG);
        store.S = S;
        store.eG = eG;
        store.G = G;
    end

    problem.hess = @hess;
    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU= reshape(YU(:)'*At*A, n, n);
        eH = 2*store.S*U + 4*sigma*(AyU*Y);
        X = Y'*Y;
        Omega = sylvester(X, X, Y'*U-Y'*U);
        T1 = U'*store.eG;
        T2 = Y'*eH;
        UY = U'*Y;
        T3 = Omega*(UY+UY');
        dOmega = sylvester(X, X, T1-T1'+T2-T2'+T3'-T3);
        H = problem.M.proj(Y, eH-U*Omega-Y*dOmega);
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
         % nY = nY/norm(nY, 'fro');
         while i <= 15 && co(nY) - cost0 > -1e-4
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              % nY = nY/norm(nY, 'fro');
              i = i + 1;
         end
    end
    
    if ~isempty(U)
        Y0 = line_search(Y0, U);
    end
    opts.verbosity = 0;
    opts.maxinner = 20;
    opts.tolgradnorm = tolgrad;
    opts.maxiter = 4;
    problem.M = quotientfactory(n, p);
    Y = trustregions(problem, Y0, opts);
end
