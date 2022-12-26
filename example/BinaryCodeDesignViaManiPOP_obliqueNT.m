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
kappa              = 6; % relaxation order

[SDP, info] = dense_sdp_relax(problem,kappa);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
Nx = K.s;


%% using mosek
% tic
% prob       = convert_sedumi2mosek(At, b,c,K);
% [~,res]    = mosekopt('minimize info',prob);
% [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
% tmosek = toc

%% using sedumi
% tic
% x = sedumi(At,b,c,K);
% tsedumi = toc
% X = reshape(x,Nx,Nx);


%% using Manopt
m = length(b);
p = 2;
A = At';
C = reshape(c, Nx, Nx);

options.maxtime = inf;
tic
[Y, fval, info] = SDP_AdptvALM_obliqueNT_subprog(A, At, b, C, c, Nx, m, p, options);
X = Y'*Y;
tmanipop = toc
fvalend = fval
msubs(f, x, X(2:d+1,1))

%disp(['Mosek:' num2str(tmosek) 's'])
%disp(['Sedumi:' num2str(tsedumi) 's'])
disp(['ManiPOP:' num2str(tmanipop) 's'])