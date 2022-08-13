%% Example: Dense SDP relaxation for binary quadratic programming (BQP)
clc; 
clear; 
close all; 

%randn('state', 11);

%% Generate binary quadratic program
d  = 4; % BQP with d variables
x  = msspoly('x',d); % symbolic decision variables using SPOTLESS
f  = (x(1)*x(2) + x(2)*x(3) + x(3)*x(4))^2 + (x(1)*x(3) + x(2)*x(4))^2 +0.0001*x(1)+0.0005*x(2)-0.0002*x(3);
%f = (x(1)*x(2) + x(2)*x(3) + x(3)*x(4))^2 + (x(1)*x(3) + x(2)*x(4))^2 ;
% objective function of the BQP
h  = x.^2 - 1; % equality constraints of the BQP (binary variables)

%% Relax BQP into an SDP
problem.vars       = x;
problem.objective  = f;
problem.equality   = h; 
kappa              = 2; % relaxation order

[SDP, info] = dense_sdp_relax(problem,kappa);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;

%% using Manopt
Nx = K.s;
p = 2;
C = reshape(c, Nx, Nx);

A = At';
m = length(b);

options.maxtime = inf;

tic
[Y, fval, info] = SDP_AdptvALM_subprog_WithOptimalCertify(A, At, b, C, c, Nx, m, p, options);
X = Y*Y';
tmanipop = toc;

%% 最优性证明，如gap为0， S1为半正定, norm(X*S1,'fro')为0
R = chol((A*At+b*b'),'lower'); %% FIX ME, maybe wrong when A*At+b*b' not PSD.
[U,SA,V] = svds(R,m);
SAD = diag(SA);
idx = find(SAD<1e-7,1)-1;
if isempty(idx)
   idx = m;
end
AATinv = U(:,1:idx)*diag((1./SAD(1:idx).^2))*U(:,1:idx)';
s = sparse(Nx*Nx,1);
while 1
    yk = AATinv*((c-s)'*At+b'*fval)';
    s = c - (yk'*A)';
    S = reshape(s, Nx, Nx);
    [V,D] = eig(S);
    S1 = V*diag(max(0,diag(D)))*V';
    S1 = S1 + S1';
    s = S1(:)/2;
    gap = b'*yk - fval;
    if abs(gap)<1e-4
        break;
    end
end

XS = trace(X*S);
minEigs = eigs(S,1,'smallestreal');
disp(['ManiPOP:' num2str(tmanipop) 's' ', gap:' num2str(gap) ', min(eig(S)):' num2str(minEigs) ', Trace(X*S):' num2str(XS)])
