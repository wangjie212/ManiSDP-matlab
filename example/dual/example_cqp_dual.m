%% Generate a random complex quadratic optimization problem
%   minimize    z'*Q*z
%   subject to  |z_i|^2 = 1,  i = 1,...,n.

% clear; 
% clc;

rng(1);
n = 5;
G = randn(n,n) + 1i*randn(n,n);
Q = (G + G')/2;                         % random Hermitian matrix

%% Generate the second-order complex SOS relaxation
[Ac, bc, cc, Kc, ~] = cqpsos(Q, [], n);

%% Convert the complex SDP to a real SDP
[A, b, c, K] = complex2real_sedumi(Ac, bc, cc, Kc);

%% Solve the realified SOS-SDP with ManiDSDP
rng(0);
clear options;
Apsd = A(:, K.f+1:end);
options.dAAt = full(sum(Apsd.^2, 2));
options.sigma0 = 1e-1;
options.sigma_min = 1e-1;
options.line_search = 1;
maxb = max(abs(b));
tic;
[~, dfval, data] = ManiDSDP_unitdiag(A, b/maxb, -c, K, options);
dtmani = toc;
dfval = dfval*maxb;
demani = max(abs([data.gap,data.pinf,data.dinf]));


%% Generate the second-order complex moment relaxation
[A, b, c, K] = bqpmom_complex_nounitdiag(Q, zeros(n,1));
[Ar, br, cr, Kr] = complex2real_sedumi(A, b, 0.5*c, K);
% [Atr, br, cr, Kr] = convertCtoR(A, b, c, K);


%% Solve the realified Moment-SDP with ManiSDP
rng(0);
options.TR_maxinner = 25;
options.line_search = 1;
tic;
[~, pfval, data] = ManiSDP_complexunitdiag(Ar', br, cr, Kr, options);
ptmani = toc;
pemani = max([data.gap, data.pinf, data.dinf]);


%% Solve with MOSEK
% prob       = convert_sedumi2mosek(A', b, -c, K);
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


% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
fprintf('ManiSDP: optimum = %.8f, eta = %.1e, time = %.2fs\n', pfval, pemani, ptmani);
fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dfval, demani, dtmani);
