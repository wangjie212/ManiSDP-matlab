%% Example: Dense SDP relaxation for binary quadratic programming (BQP)
clc; 
clear; 
close all; 

randn('state', 16);
d       = 6; % BQP with d variables
x       = sdpvar(d,1);  % 实数部分
Q       = randn(d,d); Q = Q + Q'; % a random symmetric matrix
c       = randn(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = [x.^2 == 1]; % equality constraints of the BQP (binary variables)
kappa   = 2; % relaxation order

opts = sdpsettings('solver','sparsePOP','verbose',2,'sparsePOP.reduceMomentMatSW', 1 ,... 
    'sparsePOP.relaxOrder',kappa,'sparsePOP.SDPsolverOutFile',1 , ...
    'sparsePOP.matFile', 'D:\WaveformDesignCode\hlb20220812_1.mat');
optimize(h, f, opts)

load hlb20220812_1.mat

%% Solve using MOSEK
tic
prob       = convert_sedumi2mosek(A', b,c,K);
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

[xfinal, fval] = SDP_ManiSOS_subprog(A, At, b, C, c, n, m, p,[]);

