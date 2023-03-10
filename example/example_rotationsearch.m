% spotpath   = '../../../Programs/spotless';
% addpath(genpath(spotpath));
pgdpath   = '../../STRIDE';
addpath(genpath(pgdpath));
% sdpnalpath  = '../../SDPNAL+v1.0';

%% generate random problem
rng(1);
N       = 500; % number of measurements
outrate = 0.5; % outlier rate
problem.N               = N;
problem.Covariance      = 1e-4*eye(3);
problem.v1_distribution = 'uniform';
problem.nrOutliers      = round(N*outrate);
problem.boundNoise      = true;
problem.normalize       = false;

[a, wb, R_gt, problem]   = createWahbaProblem(problem);
betasq = 1;
if problem.boundNoise; betasq = betasq*(problem.noiseBound)^2; end
if betasq == 0; betasq = 1e-2; end

SDP = QUASAR_Problem(a,wb,betasq);
[At,b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,SDP.At,SDP.C,SDP.b); 
mb = K.s;
C = full(reshape(c, mb, mb));

%% Solve using ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
tic
[~, fval, data] = ManiSDP_unittrace(At, b/(N+1), c, mb, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

% fz = [[1:length(data.fac_size)]' data.fac_size];
% residue = [[1:length(data.seta)]' log10(data.seta)];
% writematrix(fz, 'd://works/mypaper/manisdp/rs_fz_500.txt','Delimiter',' ');
% writematrix(residue, 'd://works/mypaper/manisdp/rs_residue_500.txt','Delimiter',' ');

%% Solve using STRIDE
% options.pgdStepSize     = 10; % step size, default 10
% % options.maxiterPGD      = 10; % maximum outer iterations for STRIDE, default 5-10
% options.SDPNALpath      = sdpnalpath; % provide path to SDPNAL
% options.tolADMM         = 1e-4; % tolerance for warmstart, decrease this parameter for a better warmstart (but takes more time)
% options.tolPGD          = 1e-8; % tolerance on KKT residual of the SDP
% % options.lbfgseps        = false;
% pgdopts.maxiterLBFGS    = 1000;
% pgdopts.maxiterSGS      = 300;
% options.rrOpt           = 1:3; % round the leading 3 eigenvectors to generate hypotheses
% options.rrFunName       = 'local_search_quasar'; % name of the .m file that implements the local search

% Primal initialization
% [R_gnc,info_gnc]    = GNC_Wahba(a,wb,betasq,1.4);
% q_gnc               = rotm2quat(R_gnc); q_gnc = [q_gnc(2:4),q_gnc(1)]';
% v_gnc               = kron([1;info_gnc.theta_gnc],q_gnc);
% X0                  = {v_gnc*v_gnc'};

% call STRIDE
% rng(0);
% tic
% [outPGD,X,y,S]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], options);
% by = b'*y;
% gap = abs(outPGD.pobj-by)/(abs(by)+abs(outPGD.pobj)+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% epgd = max([gap, eta, mS]);
% time_pgd = toc;

%% Solve using MOSEK
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
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval*(N+1), emani, tmani);
