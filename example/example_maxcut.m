% sdpnalpath = '../../SDPNAL+v1.0';
% addpath(genpath(sdpnalpath));

% set = ["G65", "G66", "G67", "G70", "G72", "G77", "G81"];
L = Laplacian(append('../data/Gset/', "G32", '.txt'));
C = -1/4*sparse(L);
c = C(:);

mb = size(C,1);
b = ones(mb,1);
At = sparse(mb^2, mb);
for i = 1:mb
    At((i-1)*mb+i, i) = 1;
end
K.l = 0;
K.s = mb;
blk{1,1} = 's';
blk{1,2} = mb;
nAt = sparse(mb*(mb+1)/2, mb);
for i = 1:mb
    nAt(i*(i+1)/2, i) = 1;
end

%% Solve using MOSEK
% prob       = convert_sedumi2mosek(At, b, c, K);
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

%% Solve using ManiSDP
rng(0);
clear options;
options.p0 = 40;
tic
[~, fval, data] = ManiSDP_onlyunitdiag(C, options);
emani = data.dinf;
tmani = toc;

%% Solve using SDPNAL+
% options.tol = 1e-8;
% addpath(genpath(sdpnalpath));
% rng(0);
% tic
% [objnal,X,~,y,S] = sdpnalplus(blk, {nAt}, {C}, b, [], [], [], [], [], options);
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
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
