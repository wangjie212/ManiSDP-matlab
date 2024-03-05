clear K;
[At, b, c, K] = fromsdpa('../../ManiDSDP/polyphasecode1.sdpa');
A = At';


%% Solve with ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
tic
[~, fval, data] = ManiSDP_unitdiag_multiblock1(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);