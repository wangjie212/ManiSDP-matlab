% clear; clc;
%% Generate random sparse quartic program on a sphere
rng(1);
clear I;
t = 10; % number of cliques
q = 10; % size of cliques
n = q + (q-2)*(t-1); % BQP with n variables
for i = 1:t
    I{i} = (q-2)*(i-1)+1:(q-2)*i+2;
end
sp = [];
for i = 1:t
    sp = [sp get_basis(n, 4, I{i})];
end
sp = unique(sp', 'rows');
coe = randn(size(sp, 1), 1);

%% generate moment-SDP
[At, b, c, K] = qsmom_sparse(n, I, coe);
A = At';
K.nob = 0; % This parameter indicates the first K.nob PSD cones have unit diagonals.

%% Solve with ManiSDP
rng(0);
clear options;
options.tol = 1e-4;
options.theta = 1e-4;
options.tau1 = 1e-3;
options.tau2 = 1e-2;
options.line_search = 0;
options.alpha = 0.01;
tic
[~, fval, data] = ManiSDP_multiblock(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% Solve with SDPLR
% rng(0);
% pars.feastol = 1e-4;
% nK.s = K.s;
% tic
% [x,y] = sdplr(A, b, c, nK, pars);
% vlr = c'*x;
% cy = c - At*y;
% by = b'*y;
% gap = abs(vlr-by)/(abs(by)+abs(vlr)+1);
% eta = norm(A*x - b)/(1+norm(b));
% mS = zeros(t, 1);
% ind = 1;
% for i = 1:t
%     S = reshape(cy(ind:ind+K.s(i)^2-1), K.s(i), K.s(i));
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(S, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% elr = max([eta, gap, max(mS)]);
% tlr = toc;

%% generate SOS-SDP
[A, b, c, K, dAAt] = qssos_sparse(n, I, coe);
At = A';
K.nob = 0; % This parameter indicates the first K.nob PSD cones have unit diagonals.

%% Solve with ManiDSDP
% rng(0);
% clear options;
% options.dAAt = dAAt;
% options.tol = 1e-4;
% options.theta = 1e-3;
% options.TR_maxinner = 30;
% options.tau1 = 1e-2;
% options.tau2 = 1e-1;
% options.line_search = 1;
% maxb = max(abs(b));
% tic
% [~, dval, data] = ManiDSDP_multiblock(A, b/maxb, c, K, options);
% dval = dval*maxb;
% edmani = max([data.gap, data.pinf, data.dinf]);
% tdmani = toc;

%% Solve with MOSEK
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('maximize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, K);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(sum(K.s.^2)+K.f, 1);
% x(1:K.f) = X{1};
% mS = zeros(t, 1);
% ind = K.f+1;
% for i = 1:t
%     x(ind:ind+K.s(i)^2-1) = X{i+1}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(-S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(A*x - b)/(1+norm(b));
% emosek = max([eta, gap, max(mS)]);
% tmosek = toc;

%% Solve with SDPNAL+
options.tol = 1e-4;
[blk, nAt, nC, nb] = read_sedumi(A, b, -c, K);
rng(0);
tic
[objnal,X,~,y,S] = sdpnalplus(blk, nAt, nC, nb, [], [], [], [], [], options);
by = nb'*y;
gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
x = zeros(sum(K.s.^2)+K.f, 1);
x(1:K.f) = X{1};
mS = zeros(t, 1);
ind = K.f+1;
for i = 1:t
    x(ind:ind+K.s(i)^2-1) = X{i+1}(:);
    ind = ind + K.s(i)^2;
    [~, dS] = eig(S{i}, 'vector');
    mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
end
eta = norm(A*x - nb)/(1+norm(b));
enal = max([eta, gap,  max(mS)]);
tnal = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', -objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dval, edmani, tdmani);
