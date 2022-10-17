clc; 
clear; 
close all; 
spotpath   = '../../../Programs/spotless';
addpath(genpath(spotpath));
pgdpath   = '../../STRIDE';
% sdpnalpath  = '../../SDPNAL+v1.0';
addpath(genpath(pgdpath));
%% Generate random binary quadratic program
d       = 6; % BQP with d variables
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

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);
% SDP.M       = length(info.v); % upper bound on the trace of the moment matrix
% % need the following for fast computation in the local search method
% info.v      = msspoly2degcoeff(info.v);
% info.f      = msspoly2degcoeff(info.f);
% info.J      = msspoly2degcoeff(info.J);
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
% pgdopts.pgdStepSize     = 10;
% pgdopts.SDPNALpath      = sdpnalpath;
% pgdopts.tolADMM         = 1e-4;
% pgdopts.phase1          = 1;
% pgdopts.rrOpt           = 1:3;
% pgdopts.rrFunName       = 'local_search_bqp'; % see solvers/local_search_bqp.m for implementation of local search
% pgdopts.rrPar           = info; % need the original POP formulation for local search
% pgdopts.maxiterLBFGS    = 1000;
% pgdopts.maxiterSGS      = 300;
% % pgdopts.tolLBFGS        = 1e-8;
% pgdopts.tolPGD          = 1e-6;
% pgdopts.rrFunName       = 'local_mani';
% pgdopts.rrPar.At = At;
% pgdopts.rrPar.b = b;
% pgdopts.rrPar.d = d;
% pgdopts.rrPar.c = c;
% pgdopts.rrPar.x = x;
% pgdopts.rrPar.mb = mb;
% pgdopts.rrPar.v = info.v;

% [outPGD,sXopt,syopt,sSopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
% time_pgd                    = outPGD.totaltime;
% round solutions and check optimality certificate
% res = get_performance_bqp(Xopt,yopt,Sopt,SDP,info,pgdpath);

%% Solve using MOSEK
% [At,b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,tAt,tC,tb); 
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(0)',prob);
% [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res, SDP.blk);
% tmosek = toc;
% figure; bar(eig(Xopt{1}));

%% Solve using Manopt
% m = length(b);
% p = 2;
% A = At';
% C = reshape(c, mb, mb);
% options.maxtime = inf;

% tic
% [Y, fval, info] = SDP_AdptvALM_subprog(A, At, b, C, c, mb, m, p, options);
% % X = Y'*Y;
% tmanipop = toc;
 
%% Solve using fmincon
% fobj = @(y) y'*Q*y + e'*y;
% [px,cv] = fmincon(fobj,zeros(d,1),[],[],[],[],[],[],@binary)

% addpath(genpath(sdpnalpath));
% [Y, fval] = mani_admm(SDP, At, b, c, mb, 2, sb, sB, M1, M2, 5);
% tic
% nDRS(sb, sB, M1, M2, mb, fval, sol, lb, 1e-6);
% toc
% fval = DRSPOP(At, b, c, mb, sb, sB, M1, M2);
rng(0);
tic
[Y, S, y, fval] = ALMSDP(At, b, c, mb);
toc
% tic
% [X, S, y, cx] = nALMSDP(At, b, c, mb, 1e-6);
% toc

% [SDP0.blk, SDP0.At, SDP0.C, SDP0.b, SDP0.dA] = SOStoSDP_C(f, h, x, kappa);
% tic
% [X, y, S] = spadmm(SDP0.blk, SDP0.At, SDP0.C, SDP0.b, SDP0.dA);
% toc
% [M, fval] = manisos(SDP0, At, b, c, mb, sb, sB, M1, M2, 1);

% Solve using ADMM+
% options.tol = 1e-6;
% tic
% [obj,aX,~,ay,aS] = admmplus(SDP0.blk, SDP0.At, SDP0.C, SDP0.b, [], [], [], [], [], options);
% % [obj,aX,~,ay,aS] = sdpnalplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options);
% toc

% disp(['Mosek: optimum = ' num2str(obj(1)) ', time = ' num2str(tmosek) 's'])
% disp(['Stride: optimum = ' num2str(outPGD.pobj) ', time = ' num2str(time_pgd) 's'])
% disp(['ManiPOP: optimum = ' num2str(fval) ', time = ' num2str(tmanipop) 's'])
