function [X, S, y, fval] = ALMSDPNT_GPU(At, b, c, n)
C = reshape(c, n, n);
C = gpuArray(C);
c = gpuArray(c);
At = gpuArray(At);
b = gpuArray(b);
A = At';
p = 2;
sigma = gpuArray(1e-3);
gama = gpuArray(2);
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-6;
y = gpuArray.zeros(length(b),1);
Y = gpuArray;
egrad = gpuArray.zeros(p, n);
opts = struct();
opts.verbosity = 0;      % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.mininner = 5;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm
opts.maxiter = 4;
normb = 1+norm(b);

timespend = tic;
for iter = 1:MaxIter
    Y = SDP_ALM_subprog(Y);
    X = Y'*Y;
    x = X(:);
    fval = c'*x;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(A'*y, n, n);
    DfX = C - yA;
    lamda = sum(reshape(x.*DfX(:),n,n)); % lamda = diag(DfX*X);
    S = DfX - diag(lamda);
    S = (S + S')/2; % Hermitian
    [vS, dS] = eig(S);
    v = vS(:,1)';
    mineigS = abs(dS(1));
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [UY, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > 1e-3*e(1)); % r = rank(Y)
    if r == p - 1
        q = UY(end,:);
        U = q'*v;
    elseif r < p - 1
        p = r + 1;
        Y = diag(e(1:p))*V(:,1:p)';
        U = [zeros(r,n) ; v];
    else
        U = [zeros(p,n) ; v];
        Y = [Y; zeros(1,n)];
        p = p + 1;
    end
    Y = Y + 0.1*U;
    nrms = sqrt(sum(Y.^2, 1));
    Y = bsxfun(@times, Y, 1./nrms);
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mineigS,       neta,       r,    p,    sigma,       toc(timespend));
    if max(neta, mineigS) < tao
        break;
    end
    if iter == 1 || neta > 0.7*eta
        % sigma = min(sigma*gama, 1);
        if sigma*gama > 1
            sigma = 1e-3;
        else
            sigma = sigma*gama;
        end
    end
    eta = neta;
end

    function  Y = SDP_ALM_subprog(Y0)
        problem.M = obliquefactoryNTrans(p, n);
        problem.costgrad = @costgrad;
        problem.hess = @hess;
        % Call your favorite solver.
        Y = trustregions(problem, Y0, opts);

        function [f, G, store] = costgrad(Y, store)
            X0 = Y'*Y;
            x0 = X0(:);
            Axb0 = At'*x0 - b - y/sigma;
            f = c'*x0 + 0.5*sigma*(Axb0'*Axb0);
            AxbA = A'*Axb0;
            yA0 = reshape(AxbA, n, n);
            S = C + sigma*yA0;
            egrad = 2*Y*S;            
            G = egrad - bsxfun(@times, Y, sum(Y.*egrad, 1)); % G = problem.M.egrad2rgrad(Y, eG);
        end

        % If you want to, you can specify the Riemannian Hessian as well.
        function [He, store] = hess(Y, Ydot, store)
            Xdot = Y'*Ydot;
            xdot = Xdot(:);
            AxbdotA = A'*(At'*xdot);
            yAdot = reshape(AxbdotA, n, n);
            H = 2*Ydot*S + 4*sigma*(Y*yAdot);
            He = H - bsxfun(@times, Y, sum(Y.*H, 1)) - bsxfun(@times, Ydot, sum(Y.*egrad, 1)); % He = problem.M.ehess2rhess(Y, eG, H, Ydot);
        end

        function M = obliquefactoryNTrans(n, m)
            M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
            M.dim = @() (n-1)*m;
            M.inner = @(x, d1, d2) d1(:)'*d2(:);
            M.norm = @(x, d) norm(d(:));
            M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(x - y).^2, 1)))));
            M.typicaldist = @() pi*sqrt(m);
            M.proj = @(X, U) projection(X, U);
            M.tangent = M.proj;
            M.retr = @retraction;
            % Retraction on the oblique manifold
            function y = retraction(x, d, t)
                if nargin < 3
                    td = d;
                else
                    td = t*d;
                end
                % y = normalize_columns(x + td);
                xtd = x + td;
                y = bsxfun(@times, xtd, 1./sqrt(sum((xtd).^2, 1)));
            end
            M.rand = @() random(n, m);
            M.randvec = @(x) randomvec(n, m, x);
            M.lincomb = @matrixlincomb;
            M.zerovec = @(x) zeros(n, m);
            M.transp = @(x1, x2, d) M.proj(x2, d);
        end
        function X = normalize_columns(X)
            nrms0 = sqrt(sum(X.^2, 1));
            X = bsxfun(@times, X, 1./nrms0);
        end
        function PXH = projection(X, H)
            inners = sum(X.*H, 1);
            PXH = H - bsxfun(@times, X, inners);
        end
        function x = random(n, m)
            x = normalize_columns(randn(n, m));
        end
        function d = randomvec(n, m, x)
            d = randn(n, m);
            d = projection(x, d);
            d = d / norm(d(:));
        end
    end
end
