%% Generate random sparse binary quadratic program
rng(1);
clear I;
t = 2; % number of cliques
n = 20 + 8*(t-1); % BQP with n variables
for i = 1:t
    I{i} = 8*(i-1)+1:8*i+12;
end
sp = [];
for i = 1:t
    temp = get_basis(n, 2, I{i});
    ind = true(size(temp, 2), 1);
    ind(sum(temp>1)> 0) = false;
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
coe = randn(size(sp, 1)-1, 1);

%% generate moment SDP
[At, b, c, K] = bqpmom_sparse(n, I, coe);
K.nob = length(K.s);

%% Solve using ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
options.line_search = 1;
tic
[~, fval, data] = ManiSDP_multiblock(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% Solve using MOSEK
% for i = 1:t
%     blk{i,1} = 's';
%     blk{i,2} = K.s(i);
% end
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(sum(K.s.^2), 1);
% mS = zeros(t, 1);
% ind = 1;
% for i = 1:t
%     x(ind:ind+K.s(i)^2-1) = X{i}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(At'*x - b)/(1+norm(b));
% emosek = max([eta, gap, max(mS)]);
% tmosek = toc;

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
% sdpnalpath  = '../../SDPNAL+v1.0';
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
