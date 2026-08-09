%% Generate a random complex quartic optimization problem
%   minimize    [z]_2'*Q*[z]_2
%   subject to  |z_i|^2 = 1,  i = 1,...,n.

% clear;
% clc;
rng(2);
n = 40;
N = (n+2)*(n+1)/2;
G = zeros(N, N);
G(1:n+1,1:n+1) = randn(n+1, n+1);
Q = (G + G')/2;

%% Generate the second-order complex SOS relaxation
[A, b, c, K] = cquarticsos(Q, n);

%% Solve the realified SOS-SDP with ManiDSDP
rng(0);
clear options;
Apsd = A(:, K.f+1:end);
options.dAAt = full(sum(Apsd.^2, 2));
options.kappa2 = 1e1;
options.line_search = 1;
maxb = max(abs(b));
tic;
[~, dfval, data] = ManiDSDP_unitdiag(A, b/maxb, c, K, options);
dtmani = toc;
dfval = dfval*maxb;
demani = max(abs([data.gap,data.pinf,data.dinf]));


%% Generate the second-order complex moment relaxation
[Amt, bm, cm, Km] = cquarticmom(n, Q);

%% Solve the realified Moment-SDP with ManiSDP
rng(0);
clear options;
options.sigma0 = 1e-2;
options.line_search = 1;
maxc = max(abs(cm));
tic;
[~, pfval, data] = ManiSDP_unitdiag(Amt, bm, cm/maxc, Km, options);
ptmani = toc;
pfval = pfval*maxc;
pemani = max([data.gap, data.pinf, data.dinf]);


%% Solve with MOSEK
% prob       = convert_sedumi2mosek(A', b, c, K);
% tic
% [~,res]    = mosekopt('maximize echo(3)', prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, K);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(K.s^2+K.f, 1);
% x(1:K.f) = X{1};
% x(K.f+1:K.f+K.s^2) = X{2}(:);
% [~, dS] = eig(-S{1}, 'vector');
% mS = max(0, -dS(1))/(1+abs(dS(end)));
% eta = norm(A*x - b)/(1+norm(b));
% emosek = max([eta, gap, mS]);
% tmosek = toc;


%% Solve with SDPNAL+
options.tol = 1e-8;
nc = sparse(eye(K.s));
nK.s = K.s;
[blk, nAt, nC, nb] = read_sedumi(A(2:end,2:end), b(2:end), nc(:), nK);
rng(0);
tic
[objnal,X,~,y,S] = sdpnalplus(blk, nAt, nC, nb, [], [], [], [], [], options);
by = nb'*y;
gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
x = X{1}(:);
[~, dS] = eig(S{1}, 'vector');
mS = max(0, -dS(1))/(1+abs(dS(end)));
eta = norm(A(2:end,2:end)*x - nb)/(1+norm(b));
enal = max([eta, gap, mS]);
tnal = toc;


% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', b(1)-objnal(1), enal, tnal);
% fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', pfval, pemani, ptmani);
% fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dfval, demani, dtmani);

% fprintf('&%0.1e&%0.0f&%0.1e&%0.0f&%0.1e&%0.2f&%0.1e&%0.2f\n', emosek, tmosek, enal, tnal, pemani, ptmani, demani, dtmani);

fprintf('&%0.1e&%0.0f&%0.1e&%0.0f&%0.1e&%0.0f\n', enal, tnal, pemani, ptmani, demani, dtmani);
