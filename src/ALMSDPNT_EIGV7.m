function [X, S, y, fval] = ALMSDPNT_EIGV7(At, b, c, n)
% Solve standard linear SDP problem via Augmented Lagrangian Method based 
% on Manifold Optimization
% 
% Author: WangJie, HuLiangbing
% Data: 2022-07-29
% 
% Change log:
%   HLB 2022-10-24:
%       problem using cost, grad and hess function, while ALMSDPNT_EIGV2
%       using costgrad and hess function. 
%
%   HLB 2022-10-26:
%       change neta to sqrt(Axb'*Axb), not divided by normb = 1+norm(b).
% 
%   HLB 2022-10-29:
%       1) lamda can be calculated by lamda = sum(Y0.*(Y0*DfX)); but suit for
%       small problems;
%       2) change "isempty(d)" to  "isempty(d) && neta<tao";
%       3) update sigma based on pinf and dinf(==trustregion.info(end).gradnorm);
%       4) opts.maxinner change from 20 to 180, thus gpp500-1.dat-s can be solved.


C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3; sigma_min=1e-3; sigma_max = 1e7;
gamma = 2; eta1 = 1-10*eps; eta2 = 1+10*eps; h4 = 1;
it_pinf = 0; it_dinf = 0;
MaxIter = 30000;
tolgrad = 1e-8;
tao = 1e-6;
egrad = zeros(p,n);
y = zeros(length(b),1);
gap = 1e3;
r = p - 1;
Y = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
% opts.subproblemsolver = @trs_gep; % Call your favorite solver.
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = 180;    % maximum Hessian calls per iteration
opts.mininner  = 5;
opts.maxiter   = 4;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm

timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    [Y0, ~, info]= trustregions(problem, Y, opts);
    X = Y0'*Y0;
    x = X(:);
    fval = c'*x;
    Axb = At'*x - b;
    pinf = sqrt(Axb'*Axb);
    y = y - sigma*Axb;
    yA = reshape(A'*y, n, n);
    DfX = C - yA;      
    lamda = sum(reshape(x.*DfX(:), n, n)); % for large problem, lamda = diag(DfX*X); 避免求其他乘积，从而加速
    % lamda = sum(Y0.*(Y0*DfX)); % without reshape for small problem 
    S = DfX - diag(lamda);
    S = (S + S')/2;
    [vS, dS] = eig(S);
    dS = diag(dS);
    % [vS, dS] = mineig(S);      % for large eigen decomposition problem
    % [vS, dS] = mineigimkl(S); % for small eigen decomposition problem
    dinf = info(end).gradnorm;
    d  = dS(dS<0);
    % dinf = min(abs(dS));
    if isempty(d) && pinf<tao  % S没有负的特征值，结束
       fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
                iter,    fval,       gap,       min(dS),       pinf,       dinf,       r,    p,    sigma,   toc(timespend));
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
        % Y0 = diag(e(1:r))*V(:,1:r)';  % r个主分量
        Y0 = V(:,1:r)'.*e(1:r);  % r个主分量
        p = r;
    end
    p = p + rmiusS;
    Y = [Y0; 0.1*v]; % U = [zeros(p,n) ; v]; Y = Y + 0.1*U;         
    Y = Y./sqrt(sum(Y.^2)); % nrms = sqrt(sum(Y.^2, 1)); Y = bsxfun(@times, Y, 1./nrms);

    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mineigS,       pinf,       dinf,       r,    p,    sigma,   toc(timespend));
    if max(pinf, mineigS) < tao
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

    if pinf/dinf < eta1
       it_pinf = it_pinf + 1; it_dinf = 0;
       if it_pinf >= h4
          sigma = max(sigma/gamma,sigma_min); it_pinf = 0;
       end
    elseif pinf/dinf > eta2
       it_dinf = it_dinf + 1; it_pinf = 0;
       if it_dinf >= h4
          sigma = min(sigma*gamma,sigma_max); it_dinf = 0;
       end
    end
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