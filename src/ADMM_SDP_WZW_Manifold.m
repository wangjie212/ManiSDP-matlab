function [X, y, S, fd, fp] = ADMM_SDP_WZW_Manifold(At, b, c, K, maxiter)

if nargin<5
   maxiter = 1e7;
end

%% Read SDP input data
m = length(b);
n = K.s;
p = 2;
C = reshape(c, n, n);
A = At';
normbplus1 = norm(b) + 1;
normcplus1 = norm(c) + 1;
AAT = A*At;
% invAAT = inv(AAT);
[R, flag] = chol(AAT);
icheck = 0;
while flag ~= 0 && icheck <= 100
   [R, flag] = chol(AAT + icheck*speye(m,m)*2.2204e-16); % eps = 2.2204e-16;
   icheck = icheck + 1;
   fprintf('try chol(A*At): %d times, %0.3e diag load factor.\n', icheck, icheck*2.2204e-16);
end

Atc = At'*c;

fprintf('================================= A D M M  for  S D P ===============================\n');
fprintf('Number of Linear Constraints: %3d\n',m);
fprintf('Number of Variables: %3d\n',n);

%% initializing variables
S = eye(n,n);
X = eye(n,n);
egrad = zeros(p,n);
S_V = eye(n,n);
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = 20;    % maximum Hessian calls per iteration
opts.mininner  = 5;
opts.maxiter   = 4;
opts.tolgradnorm = 1e-8; % tolerance on gradient norm
problem.costgrad = @costgrad;
problem.hess = @hess;
Y = [];

%% setting parameters
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf(' iter        objective         p-infeas         d-infeas       d-gap           mu         error\n');
fprintf('-----------------------------------------------------------------------------------------------\n');

rho = (1+sqrt(5))/2;%- .5; %rho=1.;   % step-length for the multiplier X
mu = 50; % Augmented Lagrangian penalty parameter
% Stagnation variables
gamma = 0.5; mu_min=1e-4; mu_max = 1e4; eta1 = 1; eta2 = 1; h4 = 10;
err=1e-6; it_pinf = 0; it_dinf = 0;

tic
%% ADMM main loop
Axb = At'*X(:)-b;
for count=1:maxiter

    % update y   
    %y = -invAAT*(mu.*Axb + At'*S(:) - Atc); 
    y = R\(R'\(Atc - mu.*Axb - At'*S(:)));

    % update S
    Vp = C - reshape(y'*A, n, n);
    V = Vp - mu.*X;
    problem.M = obliquefactoryNTrans(p, n);
    Y0 = trustregions(problem, Y, opts);
    S = Y0'*Y0;
    S_V = (S - V);      
    lamda = sum(reshape(S(:).*S_V(:), n, n)); % for large problem, lamda = diag(DfX*X); 避免求其他乘积，从而加速

    % update X
    X = (1/mu).*(S_V - diag(lamda));
    %X = (1-rho).*X + rho.*Xp;
    x = X(:);
    X = (X + X')/2;

    [vX, dX] = eig(X);      % for large eigen decomposition problem
    d  = dX(dX<0);
    rmiusX = min(length(d), 8); % 取小，防止Y增加太大
    v = vX(:,1:rmiusX)'; % 取主要的不大于8个负特征值对应的特征向量组
    mineigX = (d(1))/(1+dX(end));

    [~, D, VY] = svd(Y0);
    e = diag(D);
    r = sum(e > 1e-3*e(1)); % r = rank(Y)
    if r <= p - 1         
        % Y0 = diag(e(1:r))*V(:,1:r)';  % r个主分量
        Y0 = VY(:,1:r)'.*e(1:r);  % r个主分量
        p = r;
    end
    p = p + rmiusX;
    Y = [Y0; 0.1*v]; % U = [zeros(p,n) ; v]; Y = Y + 0.1*U;         
    Y = Y./sqrt(sum(Y.^2)); % nrms = sqrt(sum(Y.^2, 1)); Y = bsxfun(@times, Y, 1./nrms);  

    % calculate current error
    Axb = At'*X(:)-b;
    pinf = norm(Axb)/normbplus1;
    vs = Vp(:)-S(:);
    dinf = sqrt((vs'*vs))/normcplus1;
    cx = c'*x;
    by = b'*y;
    dgap = abs(by-cx)/(1+abs(by)+abs(cx));

    delta = max([pinf, dinf, dgap, abs(mineigX)]);
    if delta <= err
        break;
    end

    % update penalty parameter mu
    if pinf / dinf <= eta1
        it_pinf = it_pinf + 1;
        it_dinf = 0;
        if it_pinf >= h4
            mu = max(gamma * mu, mu_min);
            it_pinf = 0;
        end
    elseif pinf / dinf > eta2
        it_dinf = it_dinf + 1;
        it_pinf = 0;
        if it_dinf >= h4
            mu = min(mu / gamma, mu_max);
            it_dinf = 0;
        end
    end

    if mod(count,1)==0
        fprintf('Iter:%d, by:%.6f, cx:%.6f, pinf:%.1e, dinf:%.1e, dgap:%.1e, mu:%.3f, err:%.2e, mineigX:%0.2e\n',count,by,cx,pinf,dinf,dgap,mu,delta,mineigX);
    end
end
elapsedtime = toc;
fd = by;
fp = cx;

fprintf('Iter:%d, by:%.6f, cx:%.6f, pinf:%.1e, dinf:%.1e, dgap:%.1e, mu:%.3f, err:%.2e\n',count,by,cx,pinf,dinf,dgap,mu,delta);
fprintf('-------------------------------------------------------------------------------------\n');
lambda_min = min(eig(X)); min_el = min(min(X));
fprintf('Minimal eigenvalue of X: %4.2e; Minimal element of X:  %4.2e\n',lambda_min,min_el);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Total ADMM iterations: %3d; Final precision: %4.2e; CPU Time %2.2fsec\n',count,delta,elapsedtime);
fprintf('=====================================================================================\n');

    function [f, G, store] = costgrad(Y, store)
       S0 = Y'*Y;
       xv = S0(:) - V(:);
       f = 0.25*(xv'*xv);
       S_V = S0 - V;
       egrad = Y*S_V; 
       %G = problem.M.egrad2rgrad(Y, store.G);    
       G = egrad - Y.*sum(Y.*egrad);
    end

    function [He, store] = hess(Y, Ydot, store)
        Xdot = Y'*Ydot;
        Xdot = Xdot + Xdot';
        H = Ydot*S_V + Y*Xdot;
        He = H - Y.*sum(Y.*H) - Ydot.*sum(Y.*egrad);
    end

%     function M = obliquefactoryNTrans(n, m)
%         % M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
%         M.dim = @() (n-1)*m;
%         M.inner = @(x, d1, d2) d1(:)'*d2(:);
%         M.norm = @(x, d) norm(d(:));
%         % M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(x - y).^2, 1)))));
%         M.typicaldist = @() pi*sqrt(m);
%         M.proj = @(X, U) U - X.*sum(X.*U);
%         M.tangent = @(X, U) U - X.*sum(X.*U); %M.proj;
% 
%         M.retr = @retraction;
%         % Retraction on the oblique manifold
%         function y = retraction(x, d, t)            
%             xtd = x + d;
%             y = xtd./sqrt(sum(xtd.^2)); % y = normalize_columns(x + td);
%         end
% 
%         M.rand = @() random(n, m);
%         M.lincomb = @matrixlincomb;
%         M.zerovec = @(x) zeros(n, m);
%         M.transp = @(x1, x2, d) d - x2.*sum(x2.*d); %M.proj(x2, d);
%         M.randvec = @(x) randomvec(n, m, x);
%         % Uniform random sampling on the sphere.
%         function x = random(n, m)
%             % x = normalize_columns(randn(n, m));
%             x = randn(n, m);
%             x = x./sqrt(sum(x.^2, 1));
%         end
% 
%         function d = randomvec(n, m, x)
%             % @(X, U) U - X.*sum(X.*U);
%             d = randn(n, m);
%             d = d - x.*sum(x.*d);
%             d = d / norm(d(:));
%         end
%    end

end
