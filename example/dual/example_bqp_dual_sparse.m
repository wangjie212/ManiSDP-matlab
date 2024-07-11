%% Generate random sparse binary quadratic program
rng(1);
clear I;
t = 2; % number of cliques
n = 20 + 8*(t-1); % BQP with n variables
for i = 1:t
    I{i} = 8*(i-1)+1:8*i+12;
end
sp = [];
for i = 1:t
    temp = get_basis(n, 4, I{i});
    ind = true(size(temp, 2), 1);
    ind(sum(temp>1)> 0) = false;
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
coe = zeros(lsp, 1);
ind = true(lsp, 1);
ind(sum(sp)>2) = false;
ind(1) = false;
coe(ind) = randn(sum(ind), 1);

%% generate SOS SDP
[A, b, c, K, dAAt] = bqpsos_sparse(n, I, coe);
K.nob = length(K.s);

%% Solve using MOSEK
% prob       = convert_sedumi2mosek(A', b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, SDP.blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = X{1}(:);
% eta = norm(At'*x - b)/(1+norm(b));
% [~, dS] = eig(S{1}, 'vector');
% mS = abs(min(dS))/(1+dS(end));
% emosek = max([eta, gap, mS]);
% tmosek = toc;

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
maxb = max(abs(b));
tic
[~, fval, data] = ManiDSDP_multiblock(A, b/maxb, c, K, options);
fval = fval*maxb;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
