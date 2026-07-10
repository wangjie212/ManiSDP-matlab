clear; clc;
%% Generate random matrix completion problems
% 注意：SeDuMi 为内点法，每次迭代需对 (2m)×(2m) 的 Schur 补矩阵做 Cholesky 分解。
% m≈12600（p=100,q=100,50*n采样）时每次迭代需数分钟，属正常现象而非代码错误。
% 若需 SeDuMi 快速对比，可将 p/q 缩小至 30～40，或将 50*n 改为 3*n。
rng(2);
p = 300;
q = 300;
n = p + q;
m = 50*n;
k = 10;
M1 = (randn(p, k) + 1i*randn(p,k)) / sqrt(2);
M2 = (randn(k, q) + 1i*randn(k,q)) / sqrt(2);
M = M1 * M2;
Omega = randi([1 p*q], 1, m);
Omega = unique(Omega);
m = length(Omega);

C = eye(n);
c = C(:);
c = complex(c);

b = zeros(2*m, 1);
row = zeros(4*m, 1);
col = zeros(4*m, 1);
val = zeros(4*m, 1);

for i = 1:m
    j = ceil(Omega(i) / q);
    k_idx = mod(Omega(i), q);
    if k_idx == 0
        k_idx = q;
    end

    row1 = (j-1)*n + k_idx + p;
    row2 = (k_idx + p - 1)*n + j;

    b(2*i-1) = 2 * real(M(j, k_idx));
    row(4*i-3) = row1;
    row(4*i-2) = row2;
    col(4*i-3) = 2*i-1;
    col(4*i-2) = 2*i-1;
    val(4*i-3) = 1;
    val(4*i-2) = 1;

    b(2*i) = 2 * imag(M(j, k_idx));
    row(4*i-1) = row1;
    row(4*i)   = row2;
    col(4*i-1) = 2*i;
    col(4*i)   = 2*i;
    val(4*i-1) = -1i;
    val(4*i)   = 1i;
end

At = sparse(row, col, val, n^2, 2*m);
At = complex(At);

K.l = 0;
K.s = n;

%% ManiCSDP
clear options;
options.tol = 1e-8;
options.theta = 1e-2;
options.TR_maxinner = 6;
options.TR_maxiter = 8;
options.delta = 10;
options.tau = 1e-3;
options.alpha = 0.1;
options.p0 = k;
tic;
[~, fval, data] = ManiCSDP(At, b, c, K, options);
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;
fprintf('ManiCSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
%{
%% SeDuMi
fprintf('调用 SeDuMi (复数 SDP, 最大迭代 300 次)\n');
pars.maxiter = 300;
pars.fid     = 1;
K_sedumi = struct();
K_sedumi.s = n;         % n×n SDP 块大小
K_sedumi.scomplex = 1;  % 第 1 个块为复 Hermitian
tic;
[x_sed, y_sed, info_sed] = sedumi(At, b, c, K_sedumi, pars);
t_sedumi = toc;
if info_sed.numerr == 0
    fprintf('SeDuMi 成功完成，用时 %.2f 秒（迭代 %d 次）\n', t_sedumi, info_sed.iter);
else
    fprintf('SeDuMi 数值问题（numerr = %d），用时 %.2f 秒，当前残差：\n', info_sed.numerr, t_sedumi);
end
X_sed = reshape(x_sed, n, n);
obj_sed = real(c' * x_sed);
pinf = norm(At'*x_sed - b) / max(1, norm(b));
S_mat = reshape(c - At * y_sed, n, n);
dS    = eig((S_mat + S_mat') / 2);
dinf  = max(0, -min(real(dS))) / (1 + max(real(dS)));
gap   = abs(c'*x_sed - b'*y_sed) / (1 + abs(c'*x_sed) + abs(b'*y_sed));
fprintf('  最优值 = %.8f\n', obj_sed);
fprintf('  pinf = %.2e, dinf = %.2e, gap = %.2e\n', pinf, dinf, gap);
%}
%% convertCtoR + ManiSDP + MOSEK
max_time = 10000;
[At_r, b_r, c_r, K_r] = convertCtoR(At', b, c, K);
m_r = length(b_r);
n_r = K_r.s;
C_r = reshape(c_r, n_r, n_r);

fprintf('\n--- ManiSDP 正在求解 ---\n');
options_r=struct();

options_r.sigma0=1e-2;
options_r.sigma_min=1e-2;
options_r.theta=1e-2;
options_r.delta=10;
options_r.alpha=0.1;
options_r.TR_maxinner=40;
options_r.TR_maxiter=20;
options_r.tau1=1e-2;
options_r.tau2=5e-1;
options_r.line_search=0;

[x_r, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, options_r);


%fprintf('At_r 维度: %d × %d\n', size(At_r, 1), size(At_r, 2));
%fprintf('n_r^2 = %d, m_r = %d\n', K_r.s^2, length(b_r));


%{
fprintf('调用 MOSEK\n');
sol_mosek = struct();
try
    param.MSK_DPAR_OPTIMIZER_MAX_TIME = max_time - 10;
    param.MSK_IPAR_LOG = 0;
    prob = convert_sedumi2mosek(At_r, b_r, c_r, K_r);
    %fprintf('prob.a 维度: %d × %d\n', size(prob.a, 1), size(prob.a, 2));
    %fprintf('prob.bara 条目数: %d\n', length(prob.bara.val));
    %fprintf('prob.barc 条目数: %d\n', length(prob.barc.val));
    %fprintf('prob.blc 范数: %.4f\n',  norm(prob.blc));

    tic;
    [rcode, res] = mosekopt('minimize echo(0)', prob, param);
    t_mosek = toc;
    if rcode == 0   
        %disp(fieldnames(res.sol.itr));
        %fprintf('MOSEK pobj: %.8f\n', res.sol.itr.pobjval);
        %fprintf('MOSEK dobj: %.8f\n', res.sol.itr.dobjval);
        %fprintf('MOSEK solsta: %s\n', res.sol.itr.solsta);

        K_mosek.s = K_r.s;
        [X_mosek, y_mosek, S_mosek, mobj] = recover_mosek_sol_blk(res, K_mosek);
        if ~isempty(mobj)
            x_mosek = X_mosek{1}(:);
            pinf = norm(At_r'*x_mosek - b_r) / max(1, norm(b_r));
            by   = b_r'*y_mosek;
            gap  = abs(mobj(1) - by) / (1 + abs(mobj(1)) + abs(by));
            S_mat = S_mosek{1};
            dS    = eig(S_mat);
            dinf  = max(0, -min(dS)) / (1 + max(dS));
            fprintf('MOSEK 成功完成，用时 %.2f 秒\n', t_mosek);
            fprintf('  最优值 = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', ...
                    mobj(1), pinf, dinf, gap);
            sol_mosek.obj  = mobj(1);
            sol_mosek.time = t_mosek;
            sol_mosek.pinf = pinf;
            sol_mosek.dinf = dinf;
            sol_mosek.gap  = gap;
        else
            fprintf('MOSEK 返回空解。\n');
        end
    else
        fprintf('MOSEK 求解失败（响应码 %d）\n', rcode);
        if isfield(res, 'rmsg')
            fprintf('   MOSEK 消息: %s\n', res.rmsg);
        end
    end
catch ME
    fprintf('MOSEK 出错: %s\n', ME.message);
end
%}
