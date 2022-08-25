clc; 
clear; 
close all; 
pgdpath   = '../../STRIDE';
sdpnalpath  = '../../SDPNAL+v1.0';
addpath(genpath(pgdpath));
%% Generate random binary quadratic program
d       = 30; % BQP with d variables
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
mb = K.s;

%% Generate SOS data
[sA, dA, sB, sb] = SOStoSDP(f, h, x, kappa);
dA = 1./dA;
sB1 = [sA sB];
iD = sparse(1:length(dA),1:length(dA),dA);
iA = iD - iD*sB*(sparse(1:size(sB,2),1:size(sB,2),ones(size(sB,2),1))+sB'*iD*sB)^(-1)*sB'*iD;
M1 = sparse(1:size(sB1,2),1:size(sB1,2),ones(size(sB1,2),1)) - sB1'*iA*sB1;
M2 = sB1'*iA;

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
% pgdopts.tolLBFGS        = 1e-8;
pgdopts.tolPGD          = 1e-6;
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
% [cAt,cb,cc,cK] = SDPT3data_SEDUMIdata(SDP0.blk,SDP0.At,SDP0.C,SDP0.b); 
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

addpath(genpath(sdpnalpath));
[X, fval] = mani_admm(SDP, At, b, c, mb, sb, sB, M1, M2, 10);

% [SDP0.blk, SDP0.At, SDP0.C, SDP0.b, SDP0.dA] = SOStoSDP_C(f, h, x, kappa);
% tic
% [X, y, S] = spadmm(SDP0.blk, SDP0.At, SDP0.C, SDP0.b, SDP0.dA);
% toc
% [S, sfval] = manisos(SDP0, At, b, c, mb);

% flag = 0;
% tic 
% while flag == 0
% %% ALM方法参数设置
% sigma = 1e-3;
% gama = 13;
% % Y = full(msubs(basis, x, px));
% Y = [];
% % Y = [Y zeros(size(Y,1),1)];
% yk = zeros(m,1);
% 
% [V,D] = eig(aX{1});
% Y = [sqrt(D(end,end))*V(:,end) sqrt(D(end-1,end-1))*V(:,end-1)];
% yk = ay;

%% 迭代循环
% MaxIter = 20;
% for iter = 1:MaxIter
%     [Y, fval, info] = SDP_ALM_subprog(A, At, b, C, c, mb, p, sigma, yk, Y);
%     X = Y*Y';
%     z = X(:);
%     cx = z'*c;
%     Axb = (z'*At)' - b;
%     if norm(Axb) < 1e-4
%         break;
%     else
%         % disp(['Iter ' num2str(iter) ': fval = ' num2str(cx,10)]);
%         yk = yk + 2*Axb*sigma;
%         sigma = min(sigma*gama, 1e4);
%     end
% end
% fval = cx;
% disp(['ManiPOP: ' num2str(fval)])

% Solve using ADMM+
% options.tol = 1e-6;
% tic
% [obj,aX,~,ay,aS] = admmplus(SDP0.blk, SDP0.At, SDP0.C, SDP0.b, [], [], [], [], [], options);
% % [obj,aX,~,ay,aS] = sdpnalplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options);
% toc

% % disp(['Mosek: ' num2str(tmosek) 's'])
% disp(['ManiPOP: ' num2str(tmanipop) 's'])
% disp(['Stride: ' num2str(time_pgd) 's'])
% % disp(['Mosek: ' num2str(obj(1))])
% disp(['ManiPOP: ' num2str(fval)])
% disp(['Stride: ' num2str(outPGD.pobj)])

% [V,D] = eig(X);
% temp = zeros(length(sA), mb-1);
% for i = 1:length(sA)
%     for j = 1:mb-1
%         temp(i, j) = trace(sA{i}*V(:,j)*V(j,:));
%     end
% end
% sB0 = [temp sB];
% sol = sB0\sb;

% for i = 1:100
%     psol = sol;
%     psol(1:mb-1) = max(0,sol(1:mb-1));
%     lsol = 2*psol - sol;
%     lsol = lsol - sB'*(sB*sB')^(-1)*(sB*lsol-sb);
%     sol = sol + 1*(lsol - psol);
%     minEig = min(lsol(1:mb-1));
%     error = norm(psol - lsol);
%     disp(['step ' num2str(i) '  error:' num2str(error) ', minEig:' num2str(minEig)]);
% end

% psd = V*diag([sol(1:mb-1);0])*V';
% sol = [Mat2Vec(psd); sol(mb:end)];
% % disp(['Certify global optimality: ' num2str(tcert) 's']);
% if gap <= 1e-3
%     flag = 1;
% %    disp(['Global optimality certified!']);
% % else
% %    disp(['Global optimality not certified, use another initial point.']);
% %     p = p + 1;
% end
% end
% tmanipop = toc;
% disp(['Mosek: optimum = ' num2str(obj(1)) ', time = ' num2str(tmosek) 's'])
% disp(['Stride: optimum = ' num2str(outPGD.pobj) ', time = ' num2str(time_pgd) 's'])
% disp(['ManiPOP: optimum = ' num2str(fval) ', time = ' num2str(tmanipop) 's'])


%% Yalmip
% yx = sdpvar(d,1);
% sdpvar lower;
% p = yx'*Q*yx;
% yh = yx.^2 - 1;
% [s1,c1] = polynomial(yx,2);
% [s2,c2] = polynomial(yx,2);
% F = [sos(p-lower-[s1 s2]*yh)];
% [sol,v,sQ] = solvesos(F,-lower,[],[c1;c2;lower]);
% sQ{1} = sQ{1}([1;3;2;6;5;4],[1;3;2;6;5;4]);
% tt = zeros(15,1);
% for i = 1:15
%     tt(i) = trace(sQ{1}*sA{i});
% end
% norm(tt + sB*[value(c1);value(c2)] - sb)
% vx = monolist(yx, 2);
% clean(v{1}'*sQ{1}*v{1} + vx'*value(c1)*yh(1) + vx'*value(c2)*yh(2) + value(lower) - p,1e-6)

%% helper functions
function s = msspoly2degcoeff(f)
[~,degmat,coeff,~] = decomp(f);
s.degmat = degmat';
s.coefficient = coeff;
end
