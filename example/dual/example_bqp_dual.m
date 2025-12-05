%% Generate random binary quadratic program
rng(1)
d = 30; % BQP with d variables
Q = randn(d,d); Q = (Q + Q')/2; % a random symmetric matrix
e = randn(d,1);

%% generate moment-SDP
[At, b, c, K] = bqpmom(d, Q, e);

%% Solve with ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
tic
[~, fval, data] = ManiSDP_unitdiag(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% generate SOS-SDP
[A, b, dAAt, mb] = bqpsos(Q, e, d);

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
clear K;
K.f = 1;
K.s = mb;
c = [1; zeros(mb^2,1)];
v = zeros(size(A,1),1);
v(1) = 1;
A = [v A];
maxb = max(abs(b));
options.line_search = 1;
% options.tau2 = 1e1; % q = 80, i = 3
tic
[~, dfval, data] = ManiDSDP_unitdiag(A, b/maxb, c, K, options);
dfval = dfval*maxb;
demani = max([data.gap, data.pinf, data.dinf]);
dtmani = toc;


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


%% Solve with COPT
% problem = sedumi2copt(A, b, c, K);
% parameter.FeasTol = 1e-8;
% parameter.DualTol = 1e-8;
% parameter.RelGap = 1e-8;
% tic
% result = copt_solve(problem, parameter);
% fcopt = result.objval;
% by = b'*result.psdpi;
% gap = abs(fcopt-by)/(abs(by)+abs(fcopt)+1);
% x = zeros(K.s^2+K.f, 1);
% x(1:K.f) = result.x;
% X = zeros(K.s, K.s);
% S = zeros(K.s, K.s);
% ind = 1;
% for i = 1:K.s
%     X(i,i:K.s) = result.psdx(ind:ind+K.s-i);
%     X(i:K.s,i) = X(i,i:K.s)';
%     S(i,i:K.s) = result.psdrc(ind:ind+K.s-i);
%     S(i:K.s,i) = S(i,i:K.s)';
%     ind = ind + K.s - i + 1;
% end
% x(K.f+1:K.f+K.s^2) = X(:);
% [~, dS] = eig(-S, 'vector');
% mS = max(0, -dS(1))/(1+abs(dS(end)));
% eta = norm(A*x - b)/(1+norm(b));
% ecopt = max([eta, gap, mS]);
% tcopt = toc;


%% Solve with SDPNAL+
% options.tol = 1e-8;
% nc = sparse(eye(K.s));
% nK.s = K.s;
% [blk, nAt, nC, nb] = read_sedumi(A(2:end,2:end), b(2:end), nc(:), nK);
% rng(0);
% tic
% [objnal,X,~,y,S] = sdpnalplus(blk, nAt, nC, nb, [], [], [], [], [], options);
% by = nb'*y;
% gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
% x = X{1}(:);
% [~, dS] = eig(S{1}, 'vector');
% mS = max(0, -dS(1))/(1+abs(dS(end)));
% eta = norm(A(2:end,2:end)*x - nb)/(1+norm(b));
% enal = max([eta, gap, mS]);
% tnal = toc;


% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('CPOT: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fcopt, ecopt, tcopt);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', b(1)-objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dfval, demani, dtmani);