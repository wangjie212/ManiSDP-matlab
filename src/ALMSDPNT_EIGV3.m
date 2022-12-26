function [Y, S, y, fval, error] = ALMSDPNT_EIGV3(At, b, c, n)
% Solve standar linear SDP problem via Augmented Lagrangian Method  based 
% on Manifolds Optimization
% 
% Change log
%   HLB 2022-10-24:
%       problem using cost, grad and hess function, while ALMSDPNT_EIGV2
%       using costgrad and hess function. 

% Change log
%   WJ 2022-11-1:
%   merge the fprintf. 

C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3;
sigma_min = 1e-2; 
sigma_max = 1e7;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-8;
egrad = zeros(p, n);
y = zeros(length(b), 1);
normb = 1 + norm(b);
Y = [];
U = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
% opts.subproblemsolver = @trs_tCG_cached; % Call your favorite solver.
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.maxiter = 4;
timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    opts.tolgradnorm = tolgrad;
    [Y, ~, info] = trustregions(problem, Y, opts);
    dinf = info(end).gradnorm;
    X = Y'*Y;
    x = X(:);
    fval = c'*x;
    Axb = A*x - b;
    % neta = norm(Axb)/normb;
    pinf = sqrt(Axb'*Axb)/normb;
    neta = pinf;
    y = y - sigma*Axb;
    yA = reshape(At*y, n, n);
    eS = C - yA;
    lambda = sum(reshape(x.*eS(:), n, n));
    S = eS - diag(lambda);
    [vS, dS] = eig(S, 'vector');
    mS = abs(min(dS))/(1+dS(end));
    % sy = norm(Y*eS);
    by = b'*y + sum(lambda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > 1e-3*e(1));
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       neta,   dinf,    r,    p,    sigma,   toc(timespend));   
    error = max([neta, gap, mS]);
    if error < tao
        break;
    end
    if r <= p - 1         
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    nne = max(min(sum(dS < 0), 8), 1);
    U = [zeros(p, n); vS(:,1:nne)'];
    p = p + nne;
    Y = [Y; zeros(nne,n)]; 
    % Y = [Y; 0.1*vS(:,1:nne)'];
    % Y = Y./sqrt(sum(Y.^2));
    
%     if iter == 1 || neta > 0.7*eta
%         if sigma > 1
%             sigma = 1e-2;
%         else
%             sigma = sigma*gama;
%         end
%     end
%     eta = neta;
    
    if pinf < dinf
          sigma = max(sigma/gama, sigma_min);
    else
          sigma = min(sigma*gama, sigma_max);
    end
%    tolgrad = pinf/2;
end

    function val = co(Y)
        X = Y'*Y;
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
    end

    function nY = line_search(Y, U)
         alpha = 0.2;
         cost0 = co(Y);
         i = 1;
         nY = Y + alpha*U;
         nY = nY./sqrt(sum(nY.^2));
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              nY = nY./sqrt(sum(nY.^2));
              i = i + 1;
         end
    end

    function [f, store] = cost(Y, store)
        X0 = Y'*Y;
        x0 = X0(:);
        Axb = A*x0 - b - y/sigma;
        f = c'*x0 + 0.5*sigma*(Axb'*Axb);
    end

    function [G, store] = grad(Y, store)
        yA0 = reshape(At*Axb, n, n);
        S = C + sigma*yA0;
        egrad = 2*Y*S;
        G = egrad - Y.*sum(Y.*egrad);
    end

    % If you want to, you can specify the Riemannian Hessian as well.
    function [He, store] = hess(Y, Ydot, store)
        Xdot = Y'*Ydot;
        xdot = Xdot(:);
        AxbdotA = At*(A*xdot);
        yAdot = reshape(AxbdotA, n, n);
        H = 2*Ydot*S + 4*sigma*(Y*yAdot);
        He = H - Y.*sum(Y.*H) - Ydot.*sum(Y.*egrad);
    end

    function M = obliquefactoryNTrans(n, m)
        M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d(:));
        % M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(x - y).^2, 1)))));
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U); %M.proj;

        M.retr = @retraction;
        % Retraction on the oblique manifold
        function y = retraction(x, d)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2)); % y = normalize_columns(x + td);
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d); %M.proj(x2, d);

        % Uniform random sampling on the sphere.
        function x = random(n, m)
            % x = normalize_columns(randn(n, m));
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end
end