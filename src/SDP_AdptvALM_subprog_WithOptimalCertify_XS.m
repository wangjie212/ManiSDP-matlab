function  [xfinal, fval, info, y] = SDP_AdptvALM_subprog_WithOptimalCertify_XS(A, At, b, C, c, n, m, p, options)
    %Outer Loop Setting
    localdefaults.rho = 50;
    localdefaults.gammas = zeros(m,1); 
    localdefaults.bound = 200;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 1/3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    
    %Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-15;
    localdefaults.maxtime = inf;
    
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
    %lambdas = options.lambdas;
    y = options.gammas;
    rho = options.rho;
    oldacc = Inf;
    xCur = [];
    xPrev = xCur; 

    totaltime = tic();    
    for OuterIter = 1 : options.maxOuterIter
        %M = elliptopefactory(n, p);
        M = obliquefactory(p,n,true);
        %M = symfixedrankYYfactory(n, p);
        problem.M = M;
        problem.cost = @cost;        
        % Define the Riemannian gradient.
        problem.egrad = @egrad;
        % If you want to, you can specify the Riemannian Hessian as well.
        problem.ehess = @ehess;

        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        %inneroptions.minstepsize = options.minstepsize;

        [xCur, ~, info] = trustregions(problem, xPrev, inneroptions);         
     
        %Update Multipliers
        X = xCur*xCur';
        x = X(:);
        cx = x'*c;
        Axb = (x'*At)' - b;
        newacc = max(abs(Axb));
        y = min(options.bound, max(-options.bound,  y + rho * Axb));        

        %% It is very important !!!!!!!!!
        yA = reshape(y'*A, n, n);
        DfX = C + yA;           % Lagrangian f(X)的导数
        lamda = diag(DfX*X);    % 对角幺模diag(X)对应的对偶变量
        S = DfX - diag(lamda); % S
        by = b'*y + sum(lamda); % 对偶目标值，严格讲 by = sum(lamda) + y'*b;
        meS = min(eig(S));
        XS = norm(X*S,'fro');

        if OuterIter == 1 || newacc > options.tau * oldacc
            rho = rho/options.thetarho;
        end
        oldacc = newacc;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm); 
        
        fprintf('Iter: %d, c*x: %.8e, b*y: %0.8e, minEigS: %.5e, XS: %.5e, |Ax-b|: %0.5e, gap: %.5e, tolgradnorm: %0.8e\n', OuterIter, cx, by, meS, XS, newacc, (cx-by), tolgradnorm);
        %if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        if tolgradnorm <= options.endingtolgradnorm
            break;
        end
        if toc(totaltime) > options.maxtime
            break;
        end
        
        xPrev = xCur;
    end

    xfinal = xCur;
    fval = c'*x;

    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        cx = x'*c;
        Axb = (x'*At)' - b + y/rho;
        f = cx + 0.5*rho*(Axb'*Axb);
        AxbA = Axb'*A;
        yA = reshape(AxbA, n, n);
        S = C + rho*yA;
        store.G = 2*S*Y;
        store.S = S;
    end

    function [G, store] = egrad(Y, store)
        G = store.G;
    end

    function [H, store] = ehess(Y, Ydot, store)
        S = store.S;
        H = 2*S*Ydot;
        Xdot = Y*Ydot';
        xdot = Xdot(:);
        AxbdotA = rho*xdot'*At*A;
        AxbdotAM = reshape(AxbdotA, n, n);
        H = H + 4*AxbdotAM*Y;
    end

end