clc; 
clear; 
close all; 
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
tic
prob       = convert_sedumi2mosek(At, b,c,K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
tmosek = toc

%% Solve using Manopt
m = length(b);
p = 2;
A = At';
C = reshape(c, Nx, Nx);
%% ALM方法参数设置
sigma =0.001;
gama = 13;
Y = [];
yk = zeros(m,1);

%% 迭代循环
tic
MaxIter = 1000;
for iter = 1:MaxIter
    [Y, fval, info] = SDP_ALM_subprog(A, At, b, C, c, Nx, p, sigma, yk, Y);
    X = Y*Y';
    x = X(:);
    cx = x'*c;
    Axb = (x'*At)' - b ;
    if norm(Axb) < 1e-4
        break;
    else
        disp(['Iter:' num2str(iter) ' ,fval=' num2str(cx,10)])
        yk = yk + 2*Axb*sigma;
        sigma = sigma * gama;
        sigma = min(sigma,10000);
    end
end
fvalend = cx
tmanipop = toc;


disp(['Mosek:' num2str(tmosek) 's'])
%disp(['Sedumi:' num2str(tsedumi) 's'])
disp(['ManiPOP:' num2str(tmanipop) 's'])
disp(['Mosek:' num2str(obj(1))])
disp(['ManiPOP:' num2str(cx)])
