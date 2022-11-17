% clc; 
% clear; 
% close all; 
spotpath   = '../../../Programs/spotless';
addpath(genpath(spotpath));
pgdpath   = '../../STRIDE';
sdpnalpath  = '../../SDPNAL+v1.0';
addpath(genpath(pgdpath));
rng(1);
%% Generate random binary quadratic program
d       = 10; % BQP with d variables
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
Q       = randn(d); Q = (Q + Q')/2; % a random symmetric matrix
e       = randn(d,1);
f       = x'*Q*x + x'*e; % objective function of the BQP
h       = x.^2 - 1; % equality constraints of the BQP (binary variables)
% mon = monomials(x, 0:4);
% coe = randn(length(mon),1);
% f = coe'*mon;
% h       = sum(x.^2) - 1;
g       = []; % ask the first variable to be positive

% writematrix(Q, '../data/bqp_Q_60_3.txt');
% writematrix(e, '../data/bqp_e_60_3.txt');

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);
SDP.M       = length(info.v); % upper bound on the trace of the moment matrix need the following for fast computation in the local search method
info.v      = msspoly2degcoeff(info.v);
info.f      = msspoly2degcoeff(info.f);
info.J      = msspoly2degcoeff(info.J);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
mb = K.s;

%% Generate SOS data
% [sA, dA, sB, sb] = SOStoSDP(f, h, x, kappa);
% dA = 1./dA;
% sB1 = [sA sB];
% iD = sparse(1:length(dA),1:length(dA),dA);
% iA = iD - iD*sB*(sparse(1:size(sB,2),1:size(sB,2),ones(size(sB,2),1))+sB'*iD*sB)^(-1)*sB'*iD;
% M1 = sparse(1:size(sB1,2),1:size(sB1,2),ones(size(sB1,2),1)) - sB1'*iA*sB1;
% M2 = sB1'*iA;

%% Solve using STRIDE
pgdopts.pgdStepSize     = 10;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.tolADMM         = 1e-4;
pgdopts.phase1          = 1;
pgdopts.rrOpt           = 1:3;
pgdopts.rrFunName       = 'local_search_bqp'; % see solvers/local_search_bqp.m for implementation of local search
pgdopts.rrPar           = info; % need the original POP formulation for local search
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 300;
pgdopts.tolPGD          = 1e-8;

rng(0);
tic
[outPGD,X,y,S]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
by = b'*y;
gap = abs(outPGD.pobj-by)/(abs(by)+abs(outPGD.pobj)+1);
x = X{1}(:);
eta = norm(At'*x - b)/(1+norm(b));
[~, dS] = eig(S{1}, 'vector');
mS = abs(min(dS))/(1+dS(end));
epgd = max([gap, eta, mS]);
time_pgd = toc;

%% Solve using MOSEK
% [At,b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,tAt,tC,tb); 
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, SDP.blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% emosek = max([eta, gap, mS]);
% tmosek = toc;

%% Solve using COPT
% tic
% X0 = sdpvar(mb, mb, 'hermitian', 'real');
% F = [X0 >= 0, At'*X0(:) == b];
% obj = c'*X0(:);
% opts = sdpsettings('verbose', 1, 'solver', 'copt');
% sol = optimize(F, obj, opts);
% X = value(X0);
% S = dual(F(1));
% y = dual(F(2));
% by = b'*y;
% gap = abs(value(obj)-by)/(abs(by)+abs(value(obj))+1);
% x = X(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% ecopt = max([eta, gap, mS]);
% tcopt = toc;

%% Solve using ManiSDP
% rng(0);
% tic
% [~, ~, ~, fval, emani] = ALMSDPNT_EIGV3(At, b, c, mb);
% tmani = toc;

%% Solve using SDPNAL+
options.tol = 1e-8;
addpath(genpath(sdpnalpath));
rng(0);
tic
[objnal,X,~,y,S] = sdpnalplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options);
by = b'*y;
gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
x = X{1}(:);
eta = norm(At'*x - b)/(1+norm(b));
[~, dS] = eig(S{1}, 'vector');
mS = abs(min(dS))/(1+dS(end));
enal = max([eta, gap, mS]);
tnal = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('COPT: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vcopt, emosek, ecopt, tcopt);
fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
% fprintf('Stride: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', outPGD.pobj, epgd, time_pgd);
% fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
