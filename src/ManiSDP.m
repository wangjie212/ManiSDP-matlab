% This function solves general linear SDPs.
% Min  <C, X>
% s.t. A(X) = b,
%      X >= 0.

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
opts.maxinner = 50;     % maximum Hessian calls per iteration
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
    r = sum(e > 1e-1*e(1));
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
    nne = min(sum(dS < 0), 8);
    % U = [zeros(n, p) vS(:,1:nne)];
    p = p + nne;
    % Y = [Y zeros(n, nne)];
    Y = [Y 0.2*vS(:,1:nne)];
        
%     if iter == 1 || pinf > 0.3*opinf
%         if sigma > 1e3
%             sigma = 1e2;
%         else
%             sigma = sigma*gama;
%         end
%     end
%     opinf = pinf;
    
    if pinf < gradnorm/4
          sigma = max(sigma/gama, sigma_min);
    else
          sigma = min(sigma*gama, sigma_max);
    end
%    tolgrad = pinf;
end
   
        function nY = line_search(Y, U, t)
        X = Y'*Y;
        x = X(:);
        D = U'*U;
        vd = D(:);
        YU = Y'*U + U'*Y;
        yu = YU(:);
        p1 = c'*yu;
        p2 = c'*vd;
        q0 = b - A*x;
        q1 = A*yu;
        q2 = A*vd;
        aa = sigma/2*norm(q2)^2;
        bb = sigma*q1'*q2;
        cc = p2 - (y + sigma*q0)'*q2 + sigma/2*norm(q1)^2;
        dd = p1 - (y + sigma*q0)'*q1;
        alpha_min = 0.01;
        alpha_max = 0.2;
        sol = vpasolve(4*aa*t^3 + 3*bb*t^2 + 2*cc*t + dd == 0, t, [alpha_min alpha_max]);
        sol = eval(sol);
        v1 = aa*alpha_min^4 + bb*alpha_min^3 + cc*alpha_min^2 + dd*alpha_min;
        v2 = aa*alpha_max^4 + bb*alpha_max^3 + cc*alpha_max^2 + dd*alpha_max;
        v = [v1;v2];
        if length(sol) >= 1
            v3 = aa*sol(1)^4 + bb*sol(1)^3 + cc*sol(1)^2 + dd*sol(1);
            v = [v;v3];
        end
        if length(sol) >= 2
            v4 = aa*sol(2)^4 + bb*sol(2)^3 + cc*sol(2)^2 + dd*sol(2);
            v = [v;v4];
        end
        if length(sol) >= 3
            v5 = aa*sol(3)^4 + bb*sol(3)^3 + cc*sol(3)^2 + dd*sol(3);
            v = [v;v5];
        end
        alpha = [alpha_min;alpha_max;sol];
        [~,I] = min(v);
        alpha(I)
        nY = Y + alpha(I)*U;
        end

    function val = co(Y)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
    end

%     function nY = line_search(Y, U)
%          alpha = 0.2;
%          cost0 = co(Y);
%          i = 1;
%          nY = Y + alpha*U;
%          while i <= 15 && co(nY) - cost0 > -1e-3
%               alpha = 0.8*alpha;
%               nY = Y + alpha*U;
%               i = i + 1;
%          end
%     end

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