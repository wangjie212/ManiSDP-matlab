spotpath   = '../../../Programs/spotless';
addpath(genpath(spotpath));
pgdpath   = '../../STRIDE';
addpath(genpath(pgdpath));
sdpnalpath  = '../../SDPNAL+v1.0';

%% generate random problem
% load('../../STRIDE/data/quasar_50_1.mat') 

rng(1);
N       = 50; % number of measurements
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

%{
%% solve SDP relaxation using CVX interface and MOSEK solver
%% The syntax here is almost exactly the same as the mathematical program
%% in Proposition 8 and eq. (20) of the ICCV paper
C = cost(a,b,betasq);
n = 4*N + 4;
cvx_begin % sdp
cvx_solver mosek
variable X(n,n) symmetric
minimize( trace(C * X) ) %
subject to
X == semidefinite(n);
trace(X((1:4),(1:4))) == 1;
for k=1:N
    idx = 4+blkIndices(k,4);
    X(idx,idx) == X((1:4),(1:4));
end
for k1=1:N
    for k2=k1+1:N+1
        idx1 = blkIndices(k1,4);
        idx2 = blkIndices(k2,4);
        X(idx1,idx2) == X(idx1,idx2)';
    end
end
cvx_end

% extract solution
f_sdp   = cvx_optval; % lower bound
[V,~]   = sorteig(X);
v       = V(:,1);
q       = normalize_cols( v(1:4) );
theta   = zeros(N,1);
for i = 1:N
    theta(i) = sign(q'*v(blkIndices(i+1,4)));
end
R_err   = getAngularError(quat2rot(q),R_gt);
f_est   = cost_org(a,b,betasq,q); % upper bound
subopt  = abs(f_est - f_sdp) / (1+abs(f_est)+abs(f_sdp));

fprintf('f_sdp: %3.4e, f_est: %3.4e, R_err: %3.2e[deg].\n',f_sdp,f_est,R_err);
fprintf('Relative suboptimality: %3.2e.\n',subopt);

%}

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
% 
% % Primal initialization
% % [R_gnc,info_gnc]    = GNC_Wahba(a,wb,betasq,1.4);
% % q_gnc               = rotm2quat(R_gnc); q_gnc = [q_gnc(2:4),q_gnc(1)]';
% % v_gnc               = kron([1;info_gnc.theta_gnc],q_gnc);
% % X0                  = {v_gnc*v_gnc'};
% 
% % call STRIDE
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

%% helper function
% function Q_cost = cost(v1,v2,barc2)
% P=[1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
%    0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
%    0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
%    0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
%    -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
%    0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
%    0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
%    0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
%    -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
% P=sparse(P);
% N = size(v1,2);
% Npm = 4*N + 4;
% Q_1=zeros(Npm,Npm);
% for k=1:N
%     idx = 4+blkIndices(k,4);
%     P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
%     ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
%     Q_1((1:4),idx) = Q_1((1:4),idx)-0.5*P_k+ck/2*eye(4);
%     Q_1(idx,(1:4)) = Q_1(idx,(1:4))-0.5*P_k+ck/2*eye(4);
% end
% 
% Q_2=zeros(Npm,Npm);
% for k=1:N
% %     idx = 4+blkIndices(k,4);
%     idx = blkIndices(1,4);
%     P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
%     ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) + barc2 );
%     Q_2(idx,idx) = Q_2(idx,idx) - P_k + ck*eye(4);
% end
% Q_cost=Q_1+Q_2;
% Q_cost=sparse(Q_cost);
% 
% Q_cost=Q_cost/barc2;
% end
% 
% function B = normalize_cols(A)
% mag = sum(A.^2,1).^(0.5);
% B   = A./mag;
% end
% 
% function f_est = cost_org(a,b,betasq,q)
% R         = quat2rot(q);
% residuals = sum((b - R*a).^2,1)/betasq;
% f_est     = sum( min(residuals(:),1) );
% end

%% Solve using ManiSDP
rng(0);
tic
[Y, S, y, fval, emani] = ManiSDP_unittrace(At, b/(N+1), c, mb);
tmani = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('COPT: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vcopt, ecopt, tcopt);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
% fprintf('Stride: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', outPGD.pobj, epgd, time_pgd);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval*(N+1), emani, tmani);
