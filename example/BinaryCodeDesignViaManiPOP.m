%% Example: Dense SDP relaxation for binary quadratic programming (BQP)
clc; 
clear; 
close all; 

%% Generate random binary quadratic program
d       = 4; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
f       = (x(1)*x(2) + x(2)*x(3) + x(3)*x(4))^2 + (x(1)*x(3) + x(2)*x(4))^2 +0.00001*x(1)+0.00005*x(2)-0.00002*x(3);

 % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
g       = []; % ask the first variable to be positive

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 4; % relaxation order

[SDP,info]              = dense_sdp_relax(problem,kappa);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
Nx = K.s;


%% using mosek
tic
prob       = convert_sedumi2mosek(At, b,c,K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
tmosek = toc

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


