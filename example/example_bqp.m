% clear; clc;
% addpath(genpath('..'));
% addpath(genpath('../../mosek'));
% addpath(genpath('../../SDPLR'));
% addpath(genpath('../../spotless'));
% addpath(genpath('../../STRIDE'));
% sdpnalpath  = '../../SDPNAL+v1.0';

%% Generate random binary quadratic program
rng(1);
d       = 40; % BQP with d variables
Q       = randn(d);
Q = (Q + Q')/2; % a random symmetric matrix
e       = randn(d,1);
% x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
% f       = x'*Q*x + x'*e; % objective function of the BQP
% h       = x.^2 - 1; % equality constraints of the BQP (binary variables)

%% Relax BQP into an SDP
% problem.vars            = x;
% problem.objective       = f;
% problem.equality        = h; 
% kappa                   = 2; % relaxation order
% [SDP,info]              = dense_sdp_relax_binary(problem,kappa);
% At = SDP.sedumi.At;
% b = SDP.sedumi.b;
% c = SDP.sedumi.c;
% K = SDP.sedumi.K;
% mb = K.s;

[At, b, c, mb] = bqpmom(d, Q, e);
% K.s = mb;
% C = full(reshape(c, mb, mb));

%% Solve using ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
options.TR_maxinner = 20;
options.line_search = 1;
tic
[~, fval, data] = ManiSDP_unitdiag(At, b, c, mb, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% Solve using STRIDE
% SDP.M       = length(info.v); % upper bound on the trace of the moment matrix need the following for fast computation in the local search method
% info.v      = msspoly2degcoeff(info.v);
% info.f      = msspoly2degcoeff(info.f);
% info.J      = msspoly2degcoeff(info.J);
% pgdopts.pgdStepSize     = 10;
% pgdopts.SDPNALpath      = sdpnalpath;
% pgdopts.tolADMM         = 1e-4;
% pgdopts.phase1          = 1;
% pgdopts.rrOpt           = 1:3;
% pgdopts.rrFunName       = 'local_search_bqp'; % see solvers/local_search_bqp.m for implementation of local search
% pgdopts.rrPar           = info; % need the original POP formulation for local search
% pgdopts.maxiterLBFGS    = 1000;
% pgdopts.maxiterSGS      = 300;
% pgdopts.tolPGD          = 1e-8;
% 
% rng(0);
% tic
% [outPGD,X,y,S]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
% by = b'*y;
% gap = abs(outPGD.pobj-by)/(abs(by)+abs(outPGD.pobj)+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% epgd = max([gap, eta, mS]);
% time_pgd = toc;

%% Solve using MOSEK
% [At,b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,tAt,tC,tb); 
% prob       = convert_sedumi2mosek(At, b, c, K);
% blk = cell(1,2);
% blk{1,1} = 's';
% blk{1,2} = mb;
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% emosek = max([eta, gap, mS]);
% tmosek = toc;

%% Solve using SDPLR
% rng(0);
% pars.printlevel = 1;
% pars.feastol = 1e-8;
% tic
% [x,y] = sdplr(At', b, c, K, pars);
% vlr = c'*x;
% S = C - reshape(At*y, mb, mb);
% by = b'*y;
% gap = abs(vlr-by)/(abs(by)+abs(vlr)+1);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% elr = max([eta, gap, mS]);
% tlr = toc;

%% Solve using SDPNAL+
% options.tol = 1e-8;
% addpath(genpath(sdpnalpath));
% rng(0);
% tic
% [objnal,X,~,y,S] = sdpnalplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options);
% by = b'*y;
% gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% enal = max([eta, gap, mS]);
% tnal = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
% fprintf('Stride: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', outPGD.pobj, epgd, time_pgd);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
