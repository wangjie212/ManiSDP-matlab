%% Generate the principal components analysis problem
rng(1);
d = 20;
m = 20; % number of samples
V = rand(d, m);
k = 4;
x = msspoly('x', d);
f = 0;
for i = 1:m
    xv = x'*V(:,i);
    f = f - sum(xv.^2)/k;
end
h = [x.^3 - x; sum(x.^2) - k];

%% Generate moment relaxations
% problem.vars = x;
% problem.objective = f;
% problem.equality = h;
% [SDP, info] = dense_sdp_relax(problem, 2);

%% Solve with ManiSDP
% rng(0);
% clear options;
% options.tol = 1e-8;
% options.sigma0 = 1e1;
% options.TR_maxiter = 8;
% maxc = full(max(abs(SDP.sedumi.c)));
% tic
% [~, fval, data] = ManiSDP(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c/maxc, SDP.sedumi.K, options);
% fval = -fval*maxc;
% emani = max([data.gap, data.pinf, data.dinf]);
% tmani = toc;

%% Generate SOS relaxations
[blk, At, C, b] = sosprogram(-1, [f; -1], [], h, x, 2);
[At, b, c, K] = SDPT3data_SEDUMIdata(blk, At, C, b);

%% Solve with ManiDSDP
rng(0);
clear options;
options.tol = 1e-8;
options.sigma0 = 1;
options.TR_maxiter = 8;
maxb = full(max(abs(b)));
tic
[~, dfval, data] = ManiDSDP(At', b/maxb, -c, K, options);
dfval = -dfval*maxb;
emanid = max([data.gap, data.pinf, data.dinf]);
tmanid = toc;

%% Solve with MOSEK
% prob       = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(sum(K.s.^2)+K.f, 1);
% x(1:K.f) = X{end};
% mS = zeros(length(K.s), 1);
% ind = K.f+1;
% for i = 1:length(K.s)
%     x(ind:ind+K.s(i)^2-1) = X{i}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(At'*x - b)/(1+norm(b));
% emosek = max([eta, gap, max(mS)]);
% tmosek = toc;


% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
% fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dfval, emanid, tmanid);