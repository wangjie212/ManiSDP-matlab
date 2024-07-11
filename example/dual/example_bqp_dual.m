%% Generate random binary quadratic program
rng(1)
d = 30; % BQP with d variables
Q = randn(d,d); Q = (Q + Q')/2; % a random symmetric matrix
e = randn(d,1);

%% generate SOS SDP
[A, b, dAAt, mb] = bqpsos(Q, e, d);

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
K.f = 1;
K.s = mb;
c = [1; zeros(mb^2,1)];
v = zeros(size(A,1),1);
v(1) = 1;
A = [v A];
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP_unitdiag(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);