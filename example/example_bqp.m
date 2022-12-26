clc; 
clear; 
close all; 

%randn('state',0);
%rand('state',1);
%% Generate random binary quadratic program
d       = 30; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = rand(d); Q = (Q + Q')/2; % a random symmetric matrix
c       = rand(d,1);
f       = x'*Q*x + c'*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
g       = []; % ask the first variable to be positive

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
Nx = K.s;

%% Solve using STRIDE
%options.rrFunName       = 'local_search_quasar';
%sdpnalpath      = '../SDPNAL+v1.0'; % required for ADMM+
%options.SDPNALpath      = sdpnalpath; % provide path to SDPNAL
%[out,Xopt,yopt,Sopt] = PGDSDP(SDP.blk,SDP.At,SDP.b,SDP.C,[],options)

%% Solve using MOSEK
% tic
% prob       = convert_sedumi2mosek(At, b, c, K);
% [~,res]    = mosekopt('minimize  echo(2)', prob);
% [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
% tmosek = toc
% obj
% figure; bar(eig(Xopt{1}));

%% Solve using Manopt
m = length(b);
p = 2;
A = At';
C = reshape(c, Nx, Nx);

options.maxtime = inf;
tic
[Y, fval, info] = SDP_AdptvALM_subprog_WithOptimalCertify_LineSearch(A, At, b, C, c, Nx, m, p, options);
X = Y*Y';
tmanipop = toc
