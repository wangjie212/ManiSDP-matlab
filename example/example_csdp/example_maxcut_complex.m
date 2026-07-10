clear; clc;
cd(fileparts(mfilename('fullpath')));
% addpath(genpath('..'));
% addpath(genpath('../../mosek'));

%% Generate the max-cut problem
L = Laplacian(append('../Gset/', "G55", '.txt'));
C = -1/4*sparse(L);
[At, b, c, K] = unitdiag_constraints(C);
n = K.s;
m = length(b);
%{
%% ManiCSDP
options = struct();
options.p0 = 20;
options.AL_maxiter = 20;
options.tol = 1e-6;
options.theta = 1e-3;
options.delta = 8;
options.alpha = 1;
options.TR_maxiter = 200;
options.TR_maxinner = 100;
options.line_search = 0;
tic;
[X, obj, data] = ManiCSDP_unitdiag(At, b, c, K, options);
emani = data.dinf;
tmani = toc;
fprintf('ManiCSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', obj, emani, tmani);
%}

%% SeDuMi
fprintf('调用 SeDuMi (复数 SDP, 最大迭代 300 次)\n');
K_sedumi = struct();
K_sedumi.s = n;         % n×n SDP 块大小
K_sedumi.scomplex = 1;  % 第 1 个块为复 Hermitian（而非实对称）
pars_sedumi.maxiter = 300;
pars_sedumi.fid = 0;
tic;
try
    [x_sedumi, y_sedumi, info_sedumi] = sedumi(At, b, c, K_sedumi, pars_sedumi);
    t_sedumi = toc;
    if info_sedumi.numerr == 0
        fprintf('SeDuMi 成功完成，用时 %.2f 秒 (迭代 %d 次)\n', t_sedumi, info_sedumi.iter);
    else
        fprintf('SeDuMi 数值问题 (numerr = %d)，用时 %.2f 秒，当前残差：\n', info_sedumi.numerr, t_sedumi);
    end
    X_sedumi = reshape(x_sedumi, n, n);
    obj_sedumi = real(c' * x_sedumi);
    pinf = norm(At' * x_sedumi - b) / max(1, norm(b));
    S_mat = reshape(c - At * y_sedumi, n, n);
    dS = eig((S_mat + S_mat') / 2);
    dinf = max(0, -min(real(dS))) / (1 + max(real(dS)));
    gap = abs(c' * x_sedumi - b' * y_sedumi) / (1 + abs(c' * x_sedumi) + abs(b' * y_sedumi));
    fprintf('  最优值 = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', obj_sedumi, pinf, dinf, gap);
catch ME
    fprintf('SeDuMi 出错: %s\n', ME.message);
end

%% ManiCSDP_onlyunitdiag (complex oblique manifold, direct on C, no convertCtoR)
fprintf('调用 ManiCSDP_onlyunitdiag\n');
options_cu = struct();
options_cu.p0          = 2;
options_cu.tol         = 1e-8;
options_cu.theta       = 1e-1;
options_cu.delta       = 15;
options_cu.alpha       = 0.5;
options_cu.TR_maxiter  = 40;
options_cu.TR_maxinner = 100;
options_cu.line_search = 0;
tic;
[X_cu, obj_cu, data_cu] = ManiCSDP_onlyunitdiag(C, options_cu);
t_cu = toc;
fprintf('ManiCSDP_onlyunitdiag: optimum = %0.8f, dinf = %0.1e, time = %0.2fs\n', ...
        obj_cu, data_cu.dinf, t_cu);

%% convertCtoR + ManiSDP

[At_r, b_r, c_r, K_r] = convertCtoR(At', b, c, K);
m_r = length(b_r);
n_r = K_r.s;
C_r = reshape(c_r, n_r, n_r);
fprintf('约束转换完成，开始求解\n');

%% ManiSDP_onlycomplexunitdiag (paired-oblique manifold, no At_r needed)
fprintf('调用 ManiSDP_onlycomplexunitdiag\n');
options_pob = struct();
options_pob.p0          = 2;
options_pob.tol         = 1e-8;
options_pob.theta       = 1e-1;
options_pob.delta       = 15;
options_pob.alpha       = 0.5;
options_pob.TR_maxiter  = 40;
options_pob.TR_maxinner = 100;
options_pob.line_search = 0;
tic;
[X_pob, obj_pob, data_pob] = ManiSDP_onlycomplexunitdiag(C_r, options_pob);
t_pob = toc;
fprintf('ManiSDP_onlycomplexunitdiag: optimum = %0.8f, dinf = %0.1e, time = %0.2fs\n', ...
        obj_pob, data_pob.dinf, t_pob);
%{
fprintf('调用 ManiSDP\n');
options_r = struct();
options_r.p0          = 1;
options_r.tol         = 1e-8;
options_r.theta       = 1e-3;
options_r.delta       = 15;    % tuned: 31.8s mean, std=1.0s (Group C winner; was 18)
options_r.alpha       = 0.1;
options_r.TR_maxiter  = 6;
options_r.TR_maxinner = 40;   % tuned: 37.5s mean (Group B winner; vs 41.8s baseline)
options_r.line_search = 1;
options_r.tau1        = 1e-3;
options_r.tau2        = 1.5e-2;  % tuned: 30.1s mean, std=0.5s (fine sweep winner; Phase 4)
options_r.sigma0      = 2.0;     % tuned: Phase 1 winner (sigma0/sigma_min/gama sweep)
options_r.sigma_min   = 5e-2;     % tuned: fine sweep winner (sigma_min sweep)
options_r.gama        = 2;

[x_r, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, options_r);
%}
%{
%% convertCtoR + ManiSDP_unittrace（unit-trace 流形，利用 Tr(X_real)=n 结构）
% 变量换元：X' = X_real / n_orig，则 Tr(X')=1，对应约束 b 缩小 n_orig 倍
% 求解完成后，原始目标值 = n_orig * obj_ut_scaled
fprintf('\n调用 ManiSDP_unittrace\n');
n_orig = K.s;                % 原始复 SDP 尺寸（n）
b_ut   = b_r / n_orig;       % 缩放 b：A(X')=b/n_orig <=> A(X_real)=b

%addpath(genpath('../../ManiSDP-matlab-main/ManiSDP-matlab-main'));

options_ut = struct();
options_ut.p0          = 1;
options_ut.tol         = 1e-8;
options_ut.theta       = 1e-3;
options_ut.delta       = 15;
options_ut.alpha       = 0.05;
options_ut.TR_maxiter  = 6;
options_ut.TR_maxinner = 40;
options_ut.line_search = 1;
options_ut.tau1        = 1e-5;
options_ut.tau2        = 1e-4;
options_ut.sigma0      = 1e1;
options_ut.sigma_min   = 1e2;
options_ut.gama        = 2;

[x_ut, obj_ut_scaled, data_ut] = ManiSDP_unittrace(At_r, b_ut, c_r, K_r, options_ut);
obj_ut = n_orig * obj_ut_scaled;  % 还原为原始目标值

% 验证原始约束残差（x_real = n_orig * x_ut）
x_ut_real = n_orig * x_ut;
pinf_ut = norm(At_r' * x_ut_real(:) - b_r) / (1 + norm(b_r));
fprintf('ManiSDP_unittrace: optimum = %0.8f, pinf(original) = %0.1e, dinf = %0.1e, gap = %0.1e, time = %0.2fs\n', ...
        obj_ut, pinf_ut, data_ut.dinf, data_ut.gap, data_ut.time);
%fprintf('对比 ManiSDP:       optimum = %0.8f, pinf = %0.1e, dinf = %0.1e, gap = %0.1e, time = %0.2fs\n', ...
 %       obj_r, data_r.pinf, data_r.dinf, data_r.gap, data_r.time);
%}
%{
%% MOSEK
if n >= 3000
    fprintf('MOSEK: n=%d 过大，跳过（n>=3000 时易导致进程崩溃）\n', n);
else
fprintf('调用 MOSEK\n');
max_time = 10000;
sol_mosek = struct();
try
    prob = convert_sedumi2mosek(At_r, b_r, c_r, K_r);
    param.MSK_DPAR_OPTIMIZER_MAX_TIME = max_time - 10;
    param.MSK_IPAR_LOG = 0;
    tic;
    [rcode, res] = mosekopt('minimize echo(0)', prob, param);
    t_mosek = toc;
    if rcode == 0
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
end
%}