function [Y, S, y, fval] = ALMSDPNT_EIGV3(At, b, c, n)
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
gama = 2;
MaxIter = 300;
tolgrad = 1e-6;
tao = 1e-6;
egrad = zeros(p, n);
y = zeros(length(b), 1);
normb = 1 + norm(b);
r = p - 1;
Y = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.subproblemsolver = @trs_tCG_cached; % Call your favorite solver.
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.mininner = 5;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm
opts.maxiter = 4;

timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    Y = trustregions(problem, Y, opts);
    X = Y'*Y;
    x = X(:);
    fval = c'*x;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(A'*y, n, n);
    eS = C - yA;    
    lamda = sum(reshape(x.*eS(:), n, n)); % lamda = diag(DfX*X); 避免求其他乘积，从而加速
    S = eS - diag(lamda);
    % S = (S + S')/2;
    [vS, dS] = eig(S, 'vector');
    mS = min(dS)/(1+dS(end));
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > 1e-3*e(1));
    if r <= p - 1         
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    nne = max(min(sum(dS < 0), 4), 1);
    p = p + nne;
    Y = [Y; 0.1*vS(:,1:nne)'];        
    Y = Y./sqrt(sum(Y.^2));

    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       neta,       r,    p,    sigma,   toc(timespend));
    if max(neta, abs(mS)) < tao
        break;
    end
    if iter == 1 || neta > 0.5*eta
        % sigma = min(sigma*gama, 1);
        if sigma*gama > 1
            sigma = 1e-3;
        else
            sigma = sigma*gama;
        end
    end
    eta = neta;
end

    function [f, store] = cost(Y, store)
        X0 = Y'*Y;
        x0 = X0(:);
        Axb = At'*x0 - b - y/sigma;
        f = c'*x0 + 0.5*sigma*(Axb'*Axb);
    end

    function [G, store] = grad(Y, store)
        yA0 = reshape(A'*Axb, n, n);
        S = C + sigma*yA0;
        egrad = 2*Y*S;
        % G = egrad - bsxfun(@times, Y, sum(Y.*egrad, 1));
        G = egrad - Y.*sum(Y.*egrad);
    end

    % If you want to, you can specify the Riemannian Hessian as well.
    function [He, store] = hess(Y, Ydot, store)
        Xdot = Y'*Ydot;
        xdot = Xdot(:);
        AxbdotA = A'*(At'*xdot);
        yAdot = reshape(AxbdotA, n, n);
        H = 2*Ydot*S + 4*sigma*(Y*yAdot);
        % He = H - bsxfun(@times, Y, sum(Y.*H, 1)) - bsxfun(@times, Ydot, sum(Y.*egrad, 1)); % He = problem.M.ehess2rhess(Y, eG, H, Ydot);
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