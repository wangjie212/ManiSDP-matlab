pgdpath   = '../../STRIDE';
addpath(genpath(pgdpath));
sdpnalpath  = '../../SDPNAL+v1.0';
%% construct space of n1 x n2 hankel matrices (n1 <= n2)
rng(3);
n1 = 10;
n2 = 10;
[S,k,Scell] = hankel_struct(n1,n2);
% generate random hankel matrix
u1 = randn(k,1);
U1 = applyAffineMapCell(Scell,u1);

%% Convert the nearest hankel matrix problem to SDP
SDP     = nearest_hankel_sdp(S,u1);
[At,b,c,K] = SDPT3data_SEDUMIdata(SDP.blk,SDP.At,SDP.C,SDP.b); 
n = SDP.n;
fprintf('SDP size: n = %d, m = %d.\n\n\n',n,SDP.m);

%% Optimistic initialization using SLRA
% s.m             = n1;
% s.n             = n2;
% tol             = 0;
% [zhtls,Hu]      = htls(s.m-1,u1,'linear'); % good initialization
% opts.Rini       = zhtls';
% opts.maxiter    = 5000;
% opts.tol        = tol; 
% opts.epsrel     = tol;
% opts.epsabs     = tol; 
% opts.epsgrad    = tol;
% opts.method     = 'll'; % 'll': Levenberg–Marquardt, 'qb': BFGS
% [uslra, info]   = slra(u1, s, s.m-1, opts);
% zslra           = info.Rh'; zslra = zslra/norm(zslra);
% Uslra           = applyAffineMapCell(Scell,uslra);
% ztUnorm         = norm(zslra'*Uslra); % measure of rank deficientness
% fprintf('SLRA: norm(zt*U) = %3.2e.\n',ztUnorm);
% % lift to SDP solution
% xtld            = kron([uslra;1],zslra);
% X0              = {xtld * xtld'};

%% STRIDE
% pgdopts.pgdStepSize     = 10;
% pgdopts.maxiterSGS      = 300;
% pgdopts.maxiterLBFGS    = 1000;
% pgdopts.SDPNALpath      = sdpnalpath;
% pgdopts.tolADMM         = 5e-5;
% pgdopts.phase1          = 1;
% pgdopts.rrOpt           = 1:3;
% pgdopts.lbfgsmemory     = 10;
% pgdopts.tolPGD          = 1e-8;
% pgdopts.rrFunName       = 'rr_hankel';
% rrPar.m = n1; rrPar.n = n2; rrPar.k = k; rrPar.theta = u1; rrPar.S = S;
% pgdopts.rrPar           = rrPar;
% 
% rng(0);
% tic
% [outPGD,X,y,S] = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
% by = b'*y;
% gap = abs(outPGD.pobj-by)/(abs(by)+abs(outPGD.pobj)+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% epgd = max([gap, eta, mS]);
% time_pgd = toc;

% [z,u,f_est]     = recover_solution(Xopt,SDP.C,n1,k);
% U               = applyAffineMapCell(Scell,u);
% ztUnorm         = norm(z'*U); % measure of rank deficientness
% fprintf('norm(zt*U) = %3.2e, eta = %3.2e.\n',ztUnorm,eta);

%% Solve using ManiSDP
rng(0);
tic
[~, ~, ~, fval, emani] = ManiSDP(At, b, c, n);
tmani = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
% fprintf('Stride: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', outPGD.pobj, epgd, time_pgd);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);