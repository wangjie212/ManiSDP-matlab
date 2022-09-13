function [xfinal, fval, info] = SDP_ManiSOS_subprog1(A, At, b, C, c, n, m, p, Y0)
    %Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.x = zeros(n*n,1); % gap不为零
    localdefaults.bound = 200;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    
    %Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.maxtime = inf;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-10;
    
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
    %lambdas = options.lambdas;
    x = options.x;
    rho = options.rho;
    oldacc = Inf;

    ySCur = [];
    ySPrev = ySCur; 

    warning('off', 'manopt:getHessian:approx') 

    totaltime = tic();    
    for OuterIter = 1 : options.maxOuterIter
        %M = elliptopefactory(n, p);
        %M = symfixedrankYYfactory(n, p);
        yYoS.y = euclideanfactory(m);
        %yYoS.Y = obliquefactory(p,n,true);
        yYoS.Y = symfixedrankYYfactory(n, p);
        %manifold = obliquefactory(n,p);
        problem.M = productmanifold(yYoS);
        problem.cost = @cost;        
        % Define the Riemannian gradient.
        problem.egrad = @egrad;
        % If you want to, you can specify the Riemannian Hessian as well.
        %problem.ehess = @ehess;

        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        %inneroptions.minstepsize = options.minstepsize;

        [ySCur, fval, info] = trustregions(problem, ySPrev, inneroptions);

        yCur = ySCur.y;
        YCur = ySCur.Y;
        SCur = YCur*YCur';
        sCur = SCur(:);
        ytACur = (yCur'*A)';
        scyA = sCur - c + ytACur;
        gap = c'*x - b'*yCur;
        %Update Multipliers
        newacc = max(abs(scyA));

        if OuterIter == 1 || newacc > options.tau * oldacc
            rho = rho/options.thetarho;
        
        end
        oldacc = newacc;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm);
       
        fprintf('Iter: %d, fval: %.16e, tolgradnorm: %.16e, gap: %d\n', OuterIter, fval, tolgradnorm, gap);
        %if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        if tolgradnorm <= options.endingtolgradnorm
            break;
        end
        if toc(totaltime) > options.maxtime
            break;
        end
        x = min(options.bound, max(-options.bound, x + rho * scyA));
        ySPrev = ySCur;
    end

    xfinal = ySCur;

    % Define the cost function to be /minimized/.
    problem.cost = @cost;
    function [f, store] = cost(yYoS, store)
       y = yYoS.y;
       Y = yYoS.Y;
       S = Y*Y';
       s = S(:);
       ytA = (y'*A)';
       v = c - ytA - x/rho;
       store.V = reshape(v,n,n);
       s_v = s-v;
%        S_V = reshape(s_v,n,n);
       f = -b'*y + 0.25*rho*((s_v)'*(s_v));
%        store.YG = rho*S_V*Y; 
%        store.yG = -b + 0.5*rho*((s_v)'*At)';
%        store.S = S;
    end
 
    % Define the Riemannian gradient.
    problem.egrad = @egrad;
    function [G, store] = egrad(yYoS, store)
       y = yYoS.y;
       Y = yYoS.Y;
       S = Y*Y';
       s = S(:);
       ytA = (y'*A)';
       v = c - ytA - x/rho;
       store.V = reshape(v,n,n);
       s_v = s-v;
       S_V = reshape(s_v,n,n);
       G.Y = rho*S_V*Y; 
       G.y = -b + 0.5*rho*((s_v)'*At)';
%         G.y = store.yG;
%         G.Y = store.YG;
    end

    % If you want to, you can specify the Riemannian Hessian as well.
     problem.ehess = @ehess;
     function [H, store] = ehess(yYoS, yYoSdot, store)
         %y = yYoS.y;
         ydot = yYoSdot.y;
         Y = yYoS.Y;
         Ydot = yYoSdot.Y;
         Sdot = Y*Ydot';
         Sdot = Sdot + Sdot';
         S = store.S;
         ytAdot = (ydot'*A)';
         vdot =  c - ytAdot - x/rho;
         Vdot = reshape(vdot,n,n);
         H.Y = rho*(-store.V *Ydot - Vdot*Y + 2*S*Ydot + Ydot*(Y'*Y));
         H.y = -b + 0.5*rho*((ydot'*A)*At)' + 0.5*rho*((Sdot(:)-c+x/rho)'*At)';
     end

end