rng(1);
p = 100;
q = 100;
n = p + q;
% m = ceil(1/6*n^2);
m = 5*n;
k = 20;
M = randn(p, k)*randn(k, q);
Omega = randi([1 p*q], 1, m);
Omega = unique(Omega);
m = length(Omega);
C = eye(n);
c = C(:);
b = zeros(m,1);
At = sparse(n^2, m);
nAt = sparse(n*(n+1)/2, m);
for i = 1:m
    j = ceil(Omega(i)/q);
    k = mod(Omega(i), q);
    if k == 0
        k = q;
    end
    b(i) = 2*M(j,k);
    At((j-1)*n+k+p,i) = 1;
    At((k+p-1)*n+j,i) = 1;
    nAt((k+p)*(k+p-1)/2+j, i) = sqrt(2);
end
K.l = 0;
K.s = n;
blk{1,1} = 's';
blk{1,2} = n;

%% Solve using MOSEK
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% emosek = max([eta, gap, mS]);
% tmosek = toc;

%% Solve using SDPLR
% rng(0);
% pars.printlevel = 1;
% pars.feastol = 1e-8;
% tic
% [x,y] = sdplr(At', b, c, K, pars);
% vlr = c'*x;
% S = C - reshape(At*y, n, n);
% by = b'*y;
% gap = abs(vlr-by)/(abs(by)+abs(vlr)+1);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% elr = max([eta, gap, mS]);
% tlr = toc;

%% Solve using SDPNAL+
options.tol = 1e-8;
addpath(genpath(sdpnalpath));
rng(0);
tic
[objnal,X,~,y,S] = sdpnalplus(blk, {nAt}, {C}, b, [], [], [], [], [], options);
by = b'*y;
gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
x = X{1}(:);
eta = norm(At'*x - b)/(1+norm(b));
[~, dS] = eig(S{1}, 'vector');
mS = abs(min(dS))/(1+dS(end));
enal = max([eta, gap, mS]);
tnal = toc;

%% Solve using ManiSDP
rng(0);
tic
[~, ~, ~, fval, emani] = ManiSDP0(At, b, c, n);
tmani = toc;

% norm(X{1}(1:p,p+1:n)-M)
% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
