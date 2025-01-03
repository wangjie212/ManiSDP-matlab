%% Generate random quartic program on a sphere
rng(1);
d = 30;
coe = randn(nchoosek(d+4, 4), 1);

%% generate SOS-SDP
[A, b, c, K, dAAt] = qssos(d, coe);

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
options.delta = 6;
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);