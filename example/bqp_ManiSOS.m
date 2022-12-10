%% Example: Dense SDP relaxation for binary quadratic programming (BQP)
clc; 
clear; 
close all; 

randn('seed', 1234);
rand('state',1234);
d       = 10; % BQP with d variables
x       = sdpvar(d,1);  % 实数部分
Q       = randn(d,d); Q = (Q + Q')/2; % a random symmetric matrix
cf      = randn(d,1);
f       = x'*Q*x + cf'*x; % objective function of the BQP
h       = [x.^2 == 1]; % equality constraints of the BQP (binary variables)
kappa   = 2; % relaxation order

opts = sdpsettings('solver','sparsePOP','verbose',2,'sparsePOP.reduceMomentMatSW', 1 ,... 
    'sparsePOP.relaxOrder',kappa,'sparsePOP.SDPsolverOutFile',1 , ...
    'sparsePOP.matFile', 'D:\WaveformDesignCode\hlb20220812_1.mat');
optimize(h, f, opts)

load hlb20220812_1.mat
filename = ['hlbPOP' num2str(d) '.dat-s'];
%writesdpa(filename,A,b,c,K)

%% Solve using MOSEK
tic
prob       = convert_sedumi2mosek(A', b, c, K);
[~,res]    = mosekopt('minimize  echo(2)',prob);
blk{1,2} = K.s;
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,blk);
tmosek = toc
fp1 = obj(1)

At = A';
m = length(b);
n = K.s;
C = reshape(c,n,n);
p = 2;
n = K.s;
sigma = 2;
A = At';
R = chol(A*At);
bAAA = (R\(R'\A))'*b;
% S = sdpvar(n,n,'hermitian','real');
% opts = sdpsettings('solver','mosek','verbose', 2);
% optimize([diag(S)==1,S>=0], bAAA'*S(:) + 2*norm(-A'*(R\(R'\(A*(S(:)-c))))+S(:)-c),opts)
% S1 =  double(S)

[X, S, y, fval] = ManiSOS(At, b, c, n);
fgap = (fval-fp1)/(1+abs(fval)+abs(fp1))
% [X, y, S]= SDP_KHS(At, b, c, K);

% S = Sopt{1};
% X = Xopt{1};
% y = yopt;
% tic
% [X, y, S, fd] = ADMM_SDP_WZW_YIter(At, b, c, K);
% toc

% tic
% [X, y, S, fd] = SDPAD(At, b, c);
% toc

% [X, y, S, fp, fd] = ADMM_SDP_WZW_LowRank(A, b, c, K);
% x = S(1:d,1);
% f = x'*Q*x + cf'*x% objective function of the BQP
% fp
% fd
% fp1
%trace(S*X)

% [St,yt] = test_SDP_optimal(X, C, A, n, m);
% 
% Y = X(:,1);
% p = 1;
% [Sl,yl] = test_SDP_optimal_LowRank(Y, p, C, A, n, m)
% 
%[xfinal, fval] = SDP_AdptvALM_subprog_WithOptimalCertify_LineSearch(A, At, b, C, c, n, m, p,[]);
% [St,yt] = test_SDP_optimal(X, C, A, n, m)
