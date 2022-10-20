function [X, S, y, fval] = ALMSDPNT_EIG(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3;
gama = 2;
MaxIter = 3000;
tolgrad = 1e-8;
tao = 1e-6;
egrad = zeros(p,n);
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
problem = [];
timespend = tic;
for iter = 1:MaxIter
    Y0 = SDP_ALM_subprog(Y);
    X = Y0'*Y0;
    x = X(:);
    fval = c'*x;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(A'*y, n, n);
    DfX = C - yA;    
    lamda = sum(reshape(x.*DfX(:), n, n)); % lamda = diag(DfX*X); 避免求其他乘积，从而加速
    S = DfX - diag(lamda);
    [vS, dS] = eig(S,'vector');
    %[vS, dS] = mineig(S);
    d = dS(dS<0);
    if isempty(d)  % S没有负的特征值，结束
       fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
                iter,    fval,        gap,      min(dS),       neta,       r,    p,    sigma,   toc(timespend));

        break;
    end
    rmiusS = min(length(d), 8); % 取小，防止Y增加太大
    v = vS(:,1:rmiusS)'; % 取主要的不大于8个负特征值对应的特征向量组
    mineigS = abs(d(1));
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [~, D, V] = svd(Y0);
    e = diag(D);
    r = sum(e > 1e-3*e(1)); % r = rank(Y)
    if r <= p - 1         
        Y0 = diag(e(1:r))*V(:,1:r)';  % r个主分量
        p = r;
    end
    p = p + rmiusS;
    Y = [Y0; 0.1*v]; % U = [zeros(p,n) ; v]; Y = Y + 0.1*U;     
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
            yA0 = reshape(A'*Axb0, n, n);
            S = C + sigma*yA0;
            egrad = 2*Y*S;            
            G = egrad - bsxfun(@times, Y, sum(Y.*egrad, 1));
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
            M.egrad2rgrad = M.proj;
            M.ehess2rhess = @ehess2rhess;
            function rhess = ehess2rhess(X, egrad, ehess, U)
                PXehess = projection(X, ehess);
                inners = sum(X.*egrad, 1);
                rhess = PXehess - bsxfun(@times, U, inners);
                %rhess = PXehess - U.*reshape(inners,1,size(X,2));  % optimized by huliangbing
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
