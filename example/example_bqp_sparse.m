%% Generate random sparse binary quadratic program
rng(1);
clear I;
t = 20; % number of cliques
q = 20; % size of cliques
n = q + (q-2)*(t-1); % BQP with n variables
for i = 1:t
    I{i} = (q-2)*(i-1)+1:(q-2)*i+2;
end
sp = [];
for i = 1:t
    temp = get_basis(n, 2, I{i});
    ind = true(size(temp, 2), 1);
    ind(sum(temp>1) > 0) = false;
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
coe = randn(size(sp, 1)-1, 1);

%% generate moment-SDP
[At, b, c, K] = bqpmom_sparse(n, I, coe);
K.nob = length(K.s); % This parameter indicates the first K.nob PSD cones have unit diagonals.

%% Solve with ManiSDP
rng(0);
clear options;
options.tol = 1e-8;
options.line_search = 1;
options.tau1 = 1; % 1e1 for t = 100, q = 20
tic
[~, fval, data] = ManiSDP_multiblock(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;

%% generate SOS-SDP
sp = [];
for i = 1:t
    temp = get_basis(n, 4, I{i});
    ind = true(size(temp, 2), 1);
    ind(sum(temp>1) > 0) = false;
    sp = [sp temp(:,ind)];
end
sp = unique(sp', 'rows');
sp = sortrows(sp)';
lsp = size(sp, 2);
coe1 = zeros(lsp, 1);
ind = true(lsp, 1);
ind(sum(sp)>2) = false;
ind(1) = false;
coe1(ind) = coe;
[A, b, c, K, dAAt] = bqpsos_sparse(n, I, coe1);
At = A';
K.nob = length(K.s); % This parameter indicates the first K.nob PSD cones have unit diagonals.

%% Solve with ManiDSDP
rng(0);
clear options;
options.dAAt = dAAt;
options.tol = 1e-8;
maxb = max(abs(b));
tic
[~, dval, data] = ManiDSDP_multiblock(A, b/maxb, c, K, options);
dval = dval*maxb;
edmani = max([data.gap, data.pinf, data.dinf]);
tdmani = toc;

%% Solve with MOSEK
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('maximize echo(3)', prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, K);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(sum(K.s.^2), 1);
% mS = zeros(t, 1);
% ind = 1;
% for i = 1:t
%     x(ind:ind+K.s(i)^2-1) = X{i+1}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(-S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(A(2:end,2:end)*x - b(2:end))/(1+norm(b));
% emosek = max([eta, gap, max(mS)]);
% tmosek = toc;


%% Solve with COPT
% problem = sedumi2copt(A, b, c, K);
% parameter.FeasTol = 1e-8;
% parameter.DualTol = 1e-8;
% parameter.RelGap = 1e-8;
% tic
% result = copt_solve(problem, parameter);
% fcopt = result.objval;
% by = b'*[result.psdpi; result.pi];
% gap = abs(fcopt-by)/(abs(by)+abs(fcopt)+1);
% x = zeros(sum(K.s.^2)+K.f, 1);
% x(1:K.f) = result.x;
% mS = zeros(t, 1);
% ind = 1;
% for k = 1:t
%     X = zeros(K.s(t), K.s(t));
%     S = zeros(K.s(t), K.s(t));
%     for i = 1:K.s(t)
%         X(i,i:K.s(t)) = result.psdx(ind:ind+K.s(t)-i);
%         X(i:K.s(t),i) = X(i,i:K.s(t))';
%         S(i,i:K.s(t)) = result.psdrc(ind:ind+K.s(t)-i);
%         S(i:K.s(t),i) = S(i,i:K.s(t))';
%         ind = ind + K.s(t) - i + 1;
%     end
%     x(K.f+sum(K.s(1:k-1).^2)+1:K.f+sum(K.s(1:k).^2)) = X(:);
%     [~, dS] = eig(-S, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(A*x - b)/(1+norm(b));
% ecopt = max([eta, gap, max(mS)]);
% tcopt = toc;

%% Solve with SDPNAL+
% options.tol = 1e-8;
% nc = [];
% for i = 1:t
%     temp = sparse(eye(K.s(i)));
%     nc = [nc; temp(:)];
% end
% nK.s = K.s;
% [blk, nAt, nC, nb] = read_sedumi(A(2:end,2:end), b(2:end), nc, nK);
% rng(0);
% tic
% [objnal,X,~,y,S] = sdpnalplus(blk, nAt, nC, nb, [], [], [], [], [], options);
% by = nb'*y;
% gap = abs(objnal(1)-by)/(abs(by)+abs(objnal(1))+1);
% x = zeros(sum(K.s.^2), 1);
% mS = zeros(t, 1);
% ind = 1;
% for i = 1:t
%     x(ind:ind+K.s(i)^2-1) = X{i}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(A(2:end,2:end)*x - nb)/(1+norm(b));
% enal = max([eta, gap,  max(mS)]);
% tnal = toc;

% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('CPOT: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fcopt, ecopt, tcopt);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', -objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dval, edmani, tdmani);
