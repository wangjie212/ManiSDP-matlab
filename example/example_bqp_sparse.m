% spotpath   = '../../../Programs/spotless';
% addpath(genpath(spotpath));
% pgdpath   = '../../STRIDE';
% sdpnalpath  = '../../SDPNAL+v1.0';
% addpath(genpath(pgdpath));

%% Generate random sparse binary quadratic program
rng(1);
clear I;
t = 4; % number of cliques
n = 10 + 8*(t-1); % BQP with n variables
for i = 1:t
    I{i} = 8*(i-1)+1:8*i+2;
end
sp = [];
for i = 1:t
    temp = get_basis(n, 4, I{i});
    ind = true(size(temp, 2), 1);
    for j = 1:size(temp, 2)
        if sum(temp(:,j) > 2) > 0 || sum(mod(temp(:,j),2)) == 0
           ind(j) = false;
        end
    end
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
coe = randn(size(sp, 1), 1);
[At, b, c, K] = bqpmom_sparse(n, I, coe);
% C = full(reshape(c, mb, mb));

%% Solve using ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
tic
[~, fval, data] = ManiSDP_unitdiag_multiblock(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

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
tmosek = toc;

%% Solve using COPT
% tic
% X0 = sdpvar(mb, mb, 'hermitian', 'real');
% F = [X0 >= 0, At'*X0(:) == b];
% obj = c'*X0(:);
% opts = sdpsettings('verbose', 1, 'solver', 'sdplr');
% sol = optimize(F, obj, opts);
% X = value(X0);
% S = dual(F(1));
% y = dual(F(2));
% by = -b'*y;
% vlr = value(obj);
% gap = abs(vlr-by)/(abs(by)+abs(vlr)+1);
% x = X(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% elr = max([eta, gap, mS]);
% tlr = toc;

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
