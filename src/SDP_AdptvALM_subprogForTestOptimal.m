function  [xfinal, fval, info, gammas] = SDP_AdptvALM_subprogForTestOptimal(A, At, b, C, c, n, m, p, options)
    %Outer Loop Setting
    localdefaults.rho = 1;
    %localdefaults.gammas = b~=0; % gap恰好接近0
    localdefaults.gammas = zeros(m,1); % gap不为零
    localdefaults.bound = 200;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    
    %Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-6;
    
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
    %lambdas = options.lambdas;
    gammas = options.gammas;
    rho = options.rho;
    oldacc = Inf;

    xCur = [];
    xPrev = xCur; 

    totaltime = tic();    
    for OuterIter = 1 : options.maxOuterIter
        M = elliptopefactory(n, p);
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

        [xCur, fval, info] = trustregions(problem, xPrev, inneroptions);         
     
        %Update Multipliers
        X = xCur*xCur';
        x = X(:);
        cx = x'*c;
        Axb = (x'*At)' - b ;
        newacc = max(abs(Axb));      

% 最优性判别不对        
%         AxbA = (gammas)'*A;
%         yA = reshape(AxbA, n, n);
%         S = C + 2 *yA;  %不正定

        gap = cx - gammas'*b;
        if OuterIter == 1 || newacc > options.tau * oldacc
            rho = rho/options.thetarho;
        end
        oldacc = newacc;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm); 
       
        fprintf('Iteration: %d ,fval: %.16e, tolgradnorm: %.16e, gap: %.16e\n', OuterIter, fval, tolgradnorm, gap);
        %if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        if tolgradnorm <= options.endingtolgradnorm
            break;
        end
        if toc(totaltime) > options.maxtime
            break;
        end
        gammas = min(options.bound, max(-options.bound, gammas + rho * Axb));
        xPrev = xCur;
    end

    xfinal = xCur;

    function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        cx = x'*c;
        Axb = (x'*At)' - b + gammas/rho;
        f = cx + rho*(Axb'*Axb);
        AxbA = Axb'*A;
        yA = reshape(AxbA, n, n);
        S = C + 2 * rho*yA;
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
        AxbdotA = xdot'*At*A;
        yAdot = reshape(AxbdotA, n, n);
        H = H + 8 * rho*(yAdot*Y);
    end
end