function  [Xopt, Sopt, yopt, fopt] = SDP_AdptvALM_subprog_WithOptimalCertify_XS(A, At, b, C, c, n, m, p, options)
    %Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.gammas = zeros(m,1); 
    localdefaults.bound = 200;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 3000;
    localdefaults.numOuterItertgn = 30;
    
    %Inner Loop Setting
    localdefaults.maxInnerIter = 150;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-9;
    localdefaults.maxtime = inf;
    
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
    y = zeros(m,1); 
    rho = options.rho;
    R = [];
    totaltime = tic();
    M = obliquefactory(p,n,true);
    %M = symfixedrankYYfactory(n, p);
    problem.M = M;
    problem.cost = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;

    inneroptions.tolgradnorm = tolgradnorm;
    inneroptions.verbosity = 0;
    inneroptions.maxiter = options.maxInnerIter;
    inneroptions.useRand = 0;

    for OuterIter = 1 : options.maxOuterIter
        inneroptions.tolgradnorm = tolgradnorm;
        [R, f0, info] = trustregions(problem, R, inneroptions); 
        %[R, f0, info] = arc(problem, xPrev, inneroptions);
        
        X = R*R';
        x = X(:);
        cx = x'*c;
        Axb = (x'*At)' - b;
        newacc = max(abs(Axb));

        %Update Multipliers
        y = min(options.bound, max(-options.bound,  y + rho * Axb));        

        %% It is very important !!!!!!!!!
        yA = reshape(y'*A, n, n);
        DfX = C + yA;
        lamda = diag(DfX*X);    % 对角幺模diag(X)对应的对偶变量
        St = DfX - diag(lamda);  % S
        by = b'*y + sum(lamda); % 对偶目标值，严格讲 by = sum(lamda) + b'*y;
        
        eigsopts.isreal = true;
        meS = eigs(St, 1, 'sr', eigsopts);
        XS = norm(X*St,'fro');
        gap = abs(cx-by);

        oldacc = newacc;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm); 
        
        fprintf('Iter: %d, c*x: %.8e, b*y: %0.8e, mineigS: %.5e, XS: %.5e, |Ax-b|: %0.5e, gap: %.5e, rho:%0.5e, tolgradnorm: %0.8e\n', OuterIter, cx, by, meS, XS, newacc, gap, rho,tolgradnorm);
        if (max([XS newacc gap]) < 1e-7)  && (meS >= -1e-4) % tolgradnorm <= options.endingtolgradnorm
            break;
        end
        if toc(totaltime) > options.maxtime
            break;
        end

        if OuterIter == 1 || newacc > options.tau * oldacc
            if rho < 1e4
                rho = rho/options.thetarho;
            end
        end

    end

    Xopt = X;
    Sopt = St;
    yopt = y;
    fopt = c'*x;

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
        G = store.G ;
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
