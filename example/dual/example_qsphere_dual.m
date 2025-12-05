%% Generate random quartic program on a sphere
rng(1);
d = 10;
coe = randn(nchoosek(d+4, 4), 1);

%% generate SOS-SDP
[A, b, c, K, dAAt] = qssos(d, coe);

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
options.theta = 1e-1;
options.tau2 = 0.5;
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
options.sigma_min = 1e1;
options.delta = 1;
options.alpha = 0.001;
options.TR_maxinner = 30;
options.TR_maxiter = 6;
maxb = max(abs(b));
tic
[~, fval1, data] = nManiDSDP(A, b/maxb, c, K, options);
fval1 = fval1*maxb;
emani1 = max([data.gap, data.pinf]);
tmani1 = toc;

%% Solve with ManiDSDP
% rng(0);
% clear options;
% options.dAAt = dAAt;
% options.tol = 1e-8;
% options.sigma0 = 1e-1;
% options.theta = 1e-2;
% options.alpha = 0.001;
% options.TR_maxinner = 30;
% options.TR_maxiter = 12;
% options.tau1 = 1e-1;
% options.tau2 = 1e1;
% maxb = max(abs(b));
% tic
% [~, fval2, data] = nnManiDSDP(A, b/maxb, c, K, options);
% fval2 = fval2*maxb;
% emani2 = max([data.gap, data.pinf, data.dinf]);
% tmani2 = toc;


fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
fprintf('nManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval1, emani1, tmani1);
% fprintf('nnManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval2, emani2, tmani2);