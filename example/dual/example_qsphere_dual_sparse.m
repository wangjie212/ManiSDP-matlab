%% Generate random sparse quartic program on a sphere
rng(1);
clear I;
t = 2;
n = 10 + 8*(t-1);
for i = 1:t
    I{i} = 8*(i-1)+1:8*i+2;
end
sp = [];
for i = 1:t
    sp = [sp get_basis(n, 4, I{i})];
end
sp = unique(sp', 'rows');
coe = randn(size(sp, 1), 1);

%% generate SOS SDP
[A, b, c, K, dAAt] = qssos_sparse(n, I, coe);
K.t = zeros(1, length(K.s));
K.nob = 0;

%% Solve using MOSEK
% prob       = convert_sedumi2mosek(A', b, c, K);
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

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
options.gama = 3;
options.alpha = 0.02;
options.sigma0 = 1e-1;
options.theta = 1e-2;
options.delta = 6;
options.tao = 0.5;
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP_multiblock(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
