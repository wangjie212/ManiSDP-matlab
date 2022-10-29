function [X, y, S, fd, fp] = ADMM_SDP_WZW(A, b, c, K, maxiter)

if nargin<5
   maxiter = 1e7;
end

%% Read SDP input data
m = length(b);
n = K.s;
C = reshape(c, n, n);
At = A';

%[U,SA] = svds(A, m);
%idx = find(SA<1e-7,1)-1;
%invAAT = U(:,1:idx)*diag((1./SA(1:idx).^2))*U(:,1:idx)';
%invAAT = U*diag((1./diag(SA).^2))*U';
%invAAT = (invAAT + invAAT');
invAAT = inv(A*At);

fprintf('================================= A D M M  for  S D P ===============================\n');
fprintf('Number of Linear Constraints: %3d\n',m);
fprintf('Number of Variables: %3d\n',n);

%% initializing variables

S = eye(n,n);
X = eye(n,n);

%% setting parameters
fprintf('-------------------------------------------------------------------------------------\n');
fprintf(' iter    p-infeas     d-infeas       d-gap           mu         error      objective\n');
fprintf('-------------------------------------------------------------------------------------\n');

  
eps = 1e-6; % eps=1e-3;  final precision required
rho = (1+sqrt(5))/2;%- .5; %rho=1.;   % step-length for the multiplier X
mu = 500; % Augmented Lagrangian penalty parameter
%mu = norm(full(A),2)*1.5059400902;
% parameters for update of mu
gamma = 0.5; mu_min=1e-4; mu_max = 1e4; eta1 = 100; eta2 = 100; h4 = 100;
count = 1;  err=1; it_pinf = 0; it_dinf = 0;


tic
%% ADMM main loop
while err>eps

    % update y
    Axb = (X(:)'*At)'-b;
    y = -invAAT*(mu.*Axb + ((S(:)-c)'*At)');

    % update S
    Vp = C - reshape(y'*A,n,n);
    V = Vp - mu.*X;
    [U, d] = eig(V,'vector');
    d = max(0,d);
    S = U*diag(d)*U';
    S = (S + S')/2;

    % update X
    Xp = (1/mu).*(S - V);
    X = (1-rho).*X + rho.*Xp;
    x = X(:);

    % calculate current error
    pinf = norm(Axb); pinfs = pinf/(1+norm(b));
    vs = Vp(:)-S(:);
    dinf = sqrt((vs'*vs)); dinfs = dinf/(1+norm(c));
    cx = c'*x;
    by = b'*y;
    dgap = abs(by - cx); dgaps = dgap/(1+abs(by)+abs(cx));

    err = max(pinfs,max(dinfs,dgaps));
    count = count + 1;
    
    % update penalty parameter mu
    if pinf+dinf>0
        if pinf/dinf < eta1
            it_pinf = it_pinf + 1; it_dinf = 0;
            if it_pinf > h4
                mu = max(gamma*mu,mu_min); it_pinf = 0;
            end
        else
            if pinf/dinf > eta2
                it_dinf = it_dinf + 1; it_pinf = 0;
                if it_dinf > h4
                    mu = min(mu/gamma,mu_max); it_dinf = 0;
                end
            end
        end
    end
    
    if mod(count,100)==0
        fprintf('Iter:%d, pinf:%.3e, dinf:%.3e, dgap:%.8f, mu:%.5f, err:%.8e, by:%.8f, cx:%.8f, XS:%0.8e\n',count,pinf,dinf,dgap,mu,err,by,cx,norm(X*S,'fro'));
    end
    
    if count > maxiter
        break
    end
end
elapsedtime = toc;
fd = by;
fp = cx;

fprintf('Iter:%d, pinf:%.8f, dinf:%.8f, dgap:%.8f, mu:%.6f, err:%.8f, fd:%.8f, fp:%.8f\n',count,pinf,dinf,dgap,mu,err,by,cx);
fprintf('-------------------------------------------------------------------------------------\n');
lambda_min = min(eig(X)); min_el = min(min(X));
fprintf('Minimal eigenvalue of X: %4.2e; Minimal element of X:  %4.2e\n',lambda_min,min_el);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Total ADMM iterations: %3d; Final precision: %4.2e; CPU Time %2.2fsec\n',count,err,elapsedtime);
fprintf('=====================================================================================\n');