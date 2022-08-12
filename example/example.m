clc; 
clear; 
close all; 
pgdpath   = '../../STRIDE';
addpath(genpath(pgdpath));
%% Generate random binary quadratic program
d       = 1; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = rand(d); Q = (Q + Q')/2; % a random symmetric matrix
% e       = rand(d,1);
f       = x'*Q*x; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
g       = []; % ask the first variable to be positive

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
% basis = [1;x;monomials(x,2:kappa)];
[SDP,info]              = dense_sdp_relax(problem,kappa);
SDP.M       = length(info.v); % upper bound on the trace of the moment matrix
% need the following for fast computation in the local search method
info.v      = msspoly2degcoeff(info.v);
info.f      = msspoly2degcoeff(info.f);
info.J      = msspoly2degcoeff(info.J);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
Nx = K.s;

%% Solve using STRIDE
% sdpnalpath  = '../../SDPNAL+v1.0';
% pgdopts.pgdStepSize     = 10;
% pgdopts.SDPNALpath      = sdpnalpath;
% pgdopts.tolADMM         = 10e-5;
% pgdopts.phase1          = 1;
% pgdopts.rrOpt           = 1:3;
% pgdopts.rrFunName       = 'local_search_bqp'; % see solvers/local_search_bqp.m for implementation of local search
% pgdopts.rrPar           = info; % need the original POP formulation for local search
% pgdopts.maxiterLBFGS    = 1000;
% pgdopts.maxiterSGS      = 300;
% pgdopts.tolLBFGS        = 1e-12;
% pgdopts.tolPGD          = 1e-8;
% 
% [outPGD,sXopt,syopt,sSopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
% time_pgd                    = outPGD.totaltime;
% round solutions and check optimality certificate
% res = get_performance_bqp(Xopt,yopt,Sopt,SDP,info,pgdpath);

%% Solve using MOSEK
% tic
% prob       = convert_sedumi2mosek(At, b,c,K);
% [~,res]    = mosekopt('minimize echo(0)',prob);
% [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
% tmosek = toc;
% figure; bar(eig(Xopt{1}));

%% Solve using Manopt
m = length(b);
p = 2;
A = At';
C = reshape(c, Nx, Nx);
options.maxtime = inf;

% tic
% [Y, fval, info, gammas] = SDP_AdptvALM_subprogForTestOptimal(A, At, b, C, c, Nx, m, p, options);
% X = Y*Y';
% tmanipop = toc;

% tic
% [Y, fval, info] = SDP_AdptvALM_obliqueNT_subprog(A, At, b, C, c, Nx, m, p, options);
% X = Y'*Y;
% tmanipop = toc;

%% Solve using fmincon
% fobj = @(y) y'*Q*y + e'*y;
% [px,cv] = fmincon(fobj,zeros(d,1),[],[],[],[],[],[],@binary)

%% ALM方法参数设置
sigma = 0.001;
gama = 13;
% Y = full(msubs(basis, x, px));
Y = [];
yk = zeros(m,1);

%% 迭代循环
tic
MaxIter = 1000;
for iter = 1:MaxIter
    [Y, fval, info] = SDP_ALM_subprog(A, At, b, C, c, Nx, p, sigma, yk, Y);
    X = Y*Y';
    z = X(:);
    cx = z'*c;
    Axb = (z'*At)' - b ;
    if norm(Axb) < 1e-4
        break;
    else
        disp(['Iter:' num2str(iter) ' ,fval=' num2str(cx,10)])
        yk = yk + 2*Axb*sigma;
        sigma = sigma * gama;
        sigma = min(sigma,10000);
    end
end
fval = cx;
tmanipop = toc;
% S = C;
% for i = 1:m
%     S = S - reshape(A(i,:), Nx, Nx)*yk(i);
% end
% meig = eigs(S, 1, 'smallestreal');
% r = rank(X, 1e-4);
% up = full(msubs(f, x, X(2:d+1,1)));

% disp(['Mosek:' num2str(tmosek) 's'])
disp(['ManiPOP:' num2str(tmanipop) 's'])
%disp(['Stride:' num2str(time_pgd) 's'])
% disp(['Mosek:' num2str(obj(1))])
disp(['ManiPOP:' num2str(fval)])
%disp(['Stride:' num2str(outPGD.pobj)])
% disp(['rank:' num2str(r)])
% disp(['min eig:' num2str(meig)])
% S = C - spdiags(sum((C*Y).*Y, 2), 0, Nx, Nx);
% eigs(S, 1, 'SR')

%% SOS
sx = sym('x', [1 d]);
sf       = sx*Q*sx.' - fval;
sh = sx.^2 - 1;
[sA, sB, sb] = SOStoSDP(sf, sh, sx, kappa);
[V,D] = eig(X);
temp = zeros(length(sA), 2);
for i = 1:length(sA)
    temp(i, :) = [trace(sA{i}*V(:,1)*V(1,:)) trace(sA{i}*V(:,2)*V(2,:))];
end
sB = [sB temp];
sol = sB\sb;

%% helper functions
function s = msspoly2degcoeff(f)
[~,degmat,coeff,~] = decomp(f);
s.degmat = degmat';
s.coefficient = coeff;
end
