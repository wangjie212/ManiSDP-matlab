function [X, S, y, fval] = ManiSOS(At, b, c, n)
A = At';
R = chol(A*At);
bAAA = (R\(R'\A))'*b;
p = 2;
r = p - 1;
sigma = 50; 
MaxIter = 30000;
tolgrad = 1e-8;
tao = 1e-6;
S = zeros(n, n);
egrad = zeros(p,n);
y = zeros(length(b),1);
normc = 1 + norm(c);
stepsize = 0;
Y = [];

problem.costgrad = @costgrad;
problem.hess = @hess;
% opts.subproblemsolver = @trs_gep; % Call your favorite solver.
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = 30;    % maximum Hessian calls per iteration
opts.mininner  = 3;
opts.maxiter   = 15;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm

timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    [Y, f0, info]= trustregions(problem, Y, opts);
    S = Y'*Y;
    s = S(:);
    s_c = s - c;
    Asc = At'*s_c;
    y = -R\(R'\Asc);
    Aty = A'*y;
    Atysc = Aty + s_c;
    Atysc2 = Atysc'*Atysc;
    dinf = sqrt(Atysc2)/normc;
    fval = b'*y;
    DfSvec = bAAA + 2*sigma*Atysc;
    DfS = reshape(DfSvec, n, n);
    lamda = sum(reshape(s.*DfS(:), n, n)); % for large problem, lamda = diag(DfX*X); 避免求其他乘积，从而加速
    % lamda = sum(Y0.*(Y0*DfX)); % without reshape for small problem 
    X = DfS - diag(lamda);
    X = (X + X')/2;
    %dinf = norm(A'*y + S(:) - c)/normc;
    [vX, dX] = eig(X);
    dX = diag(dX);
    %[vS, dS] = mineig(Sopt);      % for large eigen decomposition problem
    %[vS, dS] = mineigimkl(S); % for small eigen decomposition problem
    %dinf = info(end).gradnorm;
    cx = c'*X(:);
    mineigX = abs(dX(1))/(1+abs(dX(end)));
    gap = abs(cx-fval)/(1+abs(cx)+abs(fval));
    d  = dX(dX<0);
    if ~isempty(d)
        rmiusX = min(length(d), 1); % 取小，防止Y增加太大
        v = vX(:,1:rmiusX)'; % 取主要的不大于8个负特征值对应的特征向量组
%          [~, D, V] = svd(Y);
%          e = diag(D);
%          if min(e) < 1e-6
%              break;
%          end
%         r = sum(e > 1e-3*e(1)); % r = rank(Y)
%         if r <= p - 1
%             % Y0 = diag(e(1:r))*V(:,1:r)';  % r个主分量
%             Y = V(:,1:r)'.*e(1:r);  % r个主分量
%             p = r;
%         end        
         Y = [Y; 0.1*v]; % U = [zeros(p,n) ; v]; Y = Y + 0.1*U;
         stepsize = 0;
  %       [stepsize, Y] =  linesearch_decrease(problem, [Y ; zeros(rmiusX,n)], [zeros(p,n) ; v],f0);
        p = p + rmiusX;
         Y = Y./sqrt(sum(Y.^2)); % nrms = sqrt(sum(Y.^2, 1)); Y = bsxfun(@times, Y, 1./nrms);
    end

    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigX:%0.1e, stepsize:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mineigX,       stepsize,       dinf,       r,    p,    sigma,   toc(timespend));
    if max(mineigX) < tao
        % dinf = norm(A'*y + Sopt(:) - c)/normc
        break;
    end
    %     if iter == 1 || neta > 0.7*eta
    %         % sigma = min(sigma*gama, 1);
    %         if sigma*gama > 1
    %            sigma = 1e-3;
    %         else
    %            sigma = sigma*gama;
    %         end
    %     end
    % eta = neta;

%     if pinf/dinf < eta1
%        it_pinf = it_pinf + 1; it_dinf = 0;
%        if it_pinf >= h4
%           sigma = max(sigma/gamma,sigma_min); it_pinf = 0;
%        end
%     elseif pinf/dinf > eta2 || pinf > 0.7*pinfold
%        it_dinf = it_dinf + 1; it_pinf = 0;
%        if it_dinf >= h4
%           sigma = min(sigma*gamma,sigma_max); it_dinf = 0;
%        end
%     end
%     pinfold = pinf;
end

    function [f, G, store] = costgrad(Y, store)
        S = Y'*Y;
        s = S(:);
        s_c = s - c;
        Asc = At'*s_c;
        y = -(R\(R'\Asc));
        Aty =A'*y;
        Atysc = Aty + s_c;
        Atysc2 = Atysc'*Atysc;
        f = -b'*y + sigma*Atysc2;
        DfXvex = bAAA + 2*sigma*Atysc;
        X = reshape(DfXvex, n, n);
        egrad = 2*Y*X;
        G = egrad - Y.*sum(Y.*egrad);
    end

    % If you want to, you can specify the Riemannian Hessian as well.
    function [He, store] = hess(Y, Ydot, store)
        Sdot = Y'*Ydot;
        Sdot = Sdot + Sdot';
        sdot = Sdot(:);
        Ascdot = At'*sdot;
        ydot = -(R\(R'\Ascdot));
        Atydot = A'*ydot + sdot;
        yAdot = reshape(Atydot, n, n);
        H = 2*Ydot*X + 4*sigma*(Y*yAdot);
        He = H - Y.*sum(Y.*H) - Ydot.*sum(Y.*egrad);
    end

    function M = obliquefactoryNTrans(n, m)
        % M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d(:));
        % M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(x - y).^2, 1)))));
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U); %M.proj;

        M.retr = @retraction;
        % Retraction on the oblique manifold
        function y = retraction(x, d, t)            
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