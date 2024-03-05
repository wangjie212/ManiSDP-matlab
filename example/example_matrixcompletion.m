clear; clc;
% addpath(genpath('..'));
% addpath(genpath('../../mosek'));
% addpath(genpath('../../SDPLR'));
% addpath(genpath('../../spotless'));
% addpath(genpath('../../STRIDE'));

%% Generate random matrix completion problems
rng(1);
p = 1000;
q = 1000;
n = p + q;
m = 400*n;
k = 10;
M = randn(p, k)*randn(k, q);
Omega = randi([1 p*q], 1, m);
Omega = unique(Omega);
m = length(Omega);
C = eye(n);
c = C(:);
b = zeros(m,1);
row = zeros(2*m,1);
col = zeros(2*m,1);
val = zeros(2*m,1);
nrow = zeros(m,1);
ncol = zeros(m,1);
nval = zeros(m,1);
for i = 1:m
    j = ceil(Omega(i)/q);
    k = mod(Omega(i), q);
    if k == 0
        k = q;
    end
    b(i) = 2*M(j,k);
    row(2*i-1:2*i) = [(j-1)*n+k+p; (k+p-1)*n+j];
    col(2*i-1:2*i) = [i; i];
    val(2*i-1:2*i) = [1; 1];
    nrow(i) = (k+p)*(k+p-1)/2+j;
    ncol(i) = i;
    nval(i) = sqrt(2);
end
At = sparse(row,col,val,n^2,m);
nAt = sparse(nrow,ncol,nval,n*(n+1)/2,m);
K.l = 0;
K.s = n;
blk{1,1} = 's';
blk{1,2} = n;

%% Solve using ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
options.theta = 1e-2;
options.TR_maxinner = 6;
options.TR_maxiter = 8;
options.delta = 10;
options.tao = 1e-3;
options.alpha = 0.1;
tic
[~, fval, data] = ManiSDP(At, b, c, n, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

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
% sdpnalpath  = '../../SDPNAL+v1.0';
% options.tol = 1e-8;
% addpath(genpath(sdpnalpath));
% rng(0);
% tic
% [objnal,X,~,y,S] = sdpnalplus(blk, {nAt}, {C}, b, [], [], [], [], [], options);
% by = b'*y;
% gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% enal = max([eta, gap, mS]);
% tnal = toc;

% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
