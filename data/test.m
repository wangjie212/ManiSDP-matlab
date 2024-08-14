rng(0);
clear options;
options.tol = 1e-4;
options.gama = 2;
options.alpha = 0.1;
options.sigma0 = 1e-2;
options.TR_maxinner = 50;
options.TR_maxiter = 50;
options.theta = 1e-3;
options.delta = 4;
options.tao = 1e-2;
options.line_search = 0;
SDP_1.sedumi.K.nob = 0;
tic
[~, fval, data] = ManiSDP_multiblock(SDP_1.sedumi.At, SDP_1.sedumi.b, SDP_1.sedumi.c, SDP_1.sedumi.K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;
