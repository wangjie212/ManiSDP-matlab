function [Y, S, y, fval] = ALMSDPNT(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-6;
y = zeros(length(b),1);
Y = [];
% Call your favorite solver.
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
    lamda = diag(DfX*X);
    S = DfX - diag(lamda);
    [vS, dS] = eig(S, 'vector');
    v = vS(:,1)';
    mineigS = abs(dS(1));
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [UY, e, V] = svd(Y, 'vector');
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
        iter,    fval,        gap,       mineigS,       neta,       r,    p,    sigma,   toc(timespend));
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
        Y = trustregions(problem, Y0, opts);

        function [f, G, store] = costgrad(Y, store)
            X0 = Y'*Y;
            x0 = X0(:);
            Axb0 = At'*x0 - b - y/sigma;
            f = c'*x0 + 0.5*sigma*(Axb0'*Axb0);
            AxbA = A'*Axb0;
            yA0 = reshape(AxbA, n, n);
            S0 = C + sigma*yA0;
            store.S = S0;
            store.G = 2*Y*S0;
            G = problem.M.egrad2rgrad(Y, store.G);
        end

        % If you want to, you can specify the Riemannian Hessian as well.
        function [He, store] = hess(Y, Ydot, store)
            H = 2*Ydot*store.S;
            Xdot = Y'*Ydot;
            xdot = Xdot(:);
            AxbdotA = A'*(At'*xdot);
            yAdot = reshape(AxbdotA, n, n);
            H = H + 4*sigma*(Y*yAdot);
            He = problem.M.ehess2rhess(Y, store.G, H, Ydot);
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
            M.egrad2rgrad = M.proj;
            M.ehess2rhess = @ehess2rhess;
            function rhess = ehess2rhess(X, egrad, ehess, U)
                PXehess = projection(X, ehess);
                inners = sum(X.*egrad, 1);
                rhess = PXehess - bsxfun(@times, U, inners);
                %rhess = PXehess - U.*inners;  % optimized by huliangbing
            end

            M.exp = @exponential;
            % Exponential on the oblique manifold
            function y = exponential(x, d, t)
                if nargin < 3
                    % t = 1;
                    td = d;
                else
                    td = t*d;
                end

                nrm_td = sqrt(sum(td.^2, 1));
                y = bsxfun(@times, x, cos(nrm_td)) + ...
                    bsxfun(@times, td, sin(nrm_td) ./ nrm_td);

                % For those columns where the step is 0, replace y by x
                exclude = (nrm_td == 0);
                y(:, exclude) = x(:, exclude);
            end

            M.log = @logarithm;
            function v = logarithm(x1, x2)
                v = projection(x1, x2 - x1);
                dists = real(2*asin(.5*sqrt(sum((x1 - x2).^2, 1))));
                norms = real(sqrt(sum(v.^2, 1)));
                factors = dists./norms;
                % For very close points, dists is almost equal to norms, but
                % because they are both almost zero, the division above can return
                % NaN's. To avoid that, we force those ratios to 1.
                factors(dists <= 1e-10) = 1;
                v = bsxfun(@times, v, factors);
            end

            M.retr = @retraction;
            % Retraction on the oblique manifold
            function y = retraction(x, d, t)
                if nargin < 3
                    % t = 1;
                    td = d;
                else
                    td = t*d;
                end
                y = normalize_columns(x + td);
            end

            % Inverse retraction: see spherefactory.m for background
            M.invretr = @inverse_retraction;
            function d = inverse_retraction(x, y)
                d = bsxfun(@times, y, 1./sum(x.*y, 1)) - x;
            end

            M.hash = @(x) ['z' hashmd5(x(:))];
            M.rand = @() random(n, m);
            M.randvec = @(x) randomvec(n, m, x);
            M.lincomb = @matrixlincomb;
            M.zerovec = @(x) zeros(n, m);
            M.transp = @(x1, x2, d) M.proj(x2, d);
            M.pairmean = @pairmean;
            function y = pairmean(x1, x2)
                y = x1+x2;
                y = normalize_columns(y);
            end

            vect = @(X) X(:);
            M.vec = @(x, u_mat) vect(u_mat);
            M.mat = @(x, u_vec) reshape(u_vec, [n, m]);
            M.vecmatareisometries = @() true;
        end

        % Given a matrix X, returns the same matrix but with each column scaled so
        % that they have unit 2-norm.
        function X = normalize_columns(X)
            % This is faster than norms(X, 2, 1) for small X, and as fast for large X.
            nrms0 = sqrt(sum(X.^2, 1));
            %nrms = (sum(X.^2, 1)).^0.5;
            X = bsxfun(@times, X, 1./nrms0);
            %X = X./nrms; % optimized by huliangbing 2022.02.13
        end

        % Orthogonal projection of the ambient vector H onto the tangent space at X
        function PXH = projection(X, H)
            % Compute the inner product between each vector H(:, i) with its root
            % point X(:, i), that is, X(:, i)' * H(:, i). Returns a row vector.
            inners = sum(X.*H, 1);
            PXH = H - bsxfun(@times, X, inners);
        end

        % Uniform random sampling on the sphere.
        function x = random(n, m)
            x = normalize_columns(randn(n, m));
        end

        % Random normalized tangent vector at x.
        function d = randomvec(n, m, x)
            d = randn(n, m);
            d = projection(x, d);
            d = d / norm(d(:));
        end
    end
end