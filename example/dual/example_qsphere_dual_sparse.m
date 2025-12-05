%% Generate random sparse quartic program on a sphere
rng(1);
clear I;
t = 4;
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

%% generate SOS-SDP
[A, b, c, K, dAAt] = qssos_sparse(n, I, coe);
K.nob = 0;

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
options.gama = 2;
options.alpha = 0.01;
options.sigma0 = 1e-2;
options.theta = 1e-2;
options.delta = 6;
options.tau = 0.02;
options.line_search = 0;
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP_multiblock(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
