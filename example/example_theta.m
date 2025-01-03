rng(1);
n = 1000;
m = 10*n;
Omega = randi([1 n], m, 2);
Omega = Omega(Omega(:,1) < Omega(:,2),:);
Omega = unique(Omega, 'rows');
m = length(Omega);
C = -ones(n,n);
c = C(:);
b = zeros(m,1);
b = [b; 1];
row = zeros(2*m+n,1);
col = zeros(2*m+n,1);
val = zeros(2*m+n,1);
nrow = zeros(m+n,1);
ncol = zeros(m+n,1);
nval = zeros(m+n,1);
for i = 1:m
    row(2*i-1:2*i) = [(Omega(i,1)-1)*n+Omega(i,2); (Omega(i,2)-1)*n+Omega(i,1)];
    col(2*i-1:2*i) = [i; i];
    val(2*i-1:2*i) = [1; 1];
    % nrow(i) = Omega(i,2)*(Omega(i,2)-1)/2+Omega(i,1);
    if Omega(i,2) > Omega(i,1)
        nrow(i) = Omega(i,2)*(Omega(i,2)-1)/2+Omega(i,1);
    else
        nrow(i) = Omega(i,1)*(Omega(i,1)-1)/2+Omega(i,2);
    end
    ncol(i) = i;
    nval(i) = sqrt(2);
end
for i = 1:n
    row(2*m+i) = (i-1)*n+i;
    col(2*m+i) = m+1;
    val(2*m+i) = 1;
    nrow(m+i) = i*(i-1)/2+i;
    ncol(m+i) = m+1;
    nval(m+i) = 1;
end
At = sparse(row,col,val,n^2,m+1);
nAt = sparse(nrow,ncol,nval,n*(n+1)/2,m+1);
K.l = 0;
K.s = n;
blk{1,1} = 's';
blk{1,2} = n;

%% Solve using ManiSDP
rng(0);
clear options;
options.sigma0 = 1e-1;
options.tol = 1e-6;
options.TR_maxinner = 20;
options.TR_maxiter = 4;
options.tao = 1e-3;
options.theta = 1e-3;
options.delta = 8;
options.line_search = 0;
options.alpha = 0.001;
tic
[~, fval, data] = ManiSDP_unittrace(At, b, c, n, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% Solve using MOSEK
prob       = convert_sedumi2mosek(At, b, c, K);
tic
[~,res]    = mosekopt('minimize echo(3)',prob);
[X,y,S,mobj] = recover_mosek_sol_blk(res, K);
by = b'*y;
gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
x = X{1}(:);
eta = norm(At'*x - b)/(1+norm(b));
[~, dS] = eig(S{1}, 'vector');
mS = abs(min(dS))/(1+dS(end));
emosek = max([eta, gap, mS]);
tmosek = toc;

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
% addpath(genpath(sdpnalpath));
% rng(0);
% options.tol = 1e-8;
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

fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
