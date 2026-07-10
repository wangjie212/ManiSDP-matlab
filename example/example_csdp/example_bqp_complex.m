
clear; 
clc;

%% 1. 随机生成 BQP 问题参数
rng(1);
d = 35;
Q = (randn(d) + 1i*randn(d)) / sqrt(2);
Q = (Q + Q') / 2;
e = (randn(d, 1) + 1i*randn(d, 1)) / sqrt(2);

%% 2. 生成 BQP 矩量松弛数据
fprintf('正在生成 BQP 矩量松弛数据...\n');
% 含单位对角约束（供 SeDuMi/MOSEK，需完整约束）
[At_c, b_c, c_c, K_c]         = bqpmom_complex(Q, e);
% 无单位对角约束（供 ManiCSDP/ManiSDP_complexunitdiag，由流形隐式处理）
[At_cnd, b_cnd, c_cnd, K_cnd] = bqpmom_complex_nounitdiag(Q, e);
N_c = K_cnd.s;
%}
%% 3. ManiCSDP_unitdiag (直接求解复数 SDP)
fprintf('\n=== 情形 1: ManiCSDP_unitdiag (复数 SDP) ===\n');
options_c.tol         = 1e-8;
options_c.AL_maxiter  = 1000;
options_c.p0          = 2;
options_c.sigma0      = 1e-3;
options_c.sigma_min   = 1e-2;
options_c.sigma_max   = 1e7;
options_c.gama        = 2;
options_c.tau1        = 1;
options_c.tau2        = 1;
options_c.theta       = 1e-3;
options_c.delta       = 8;
options_c.alpha       = 0.1;
options_c.TR_maxinner = 25; 
options_c.TR_maxiter  = 4;
options_c.line_search = 0;

[X_c, fval_c, data_c] = ManiCSDP_unitdiag(At_cnd', b_cnd, c_cnd, K_cnd, options_c);
fprintf('ManiCSDP_unitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_c, max([data_c.gap, data_c.pinf, data_c.dinf]), data_c.time);
%{
%% 4. SeDuMi (直接求解复数 SDP，含完整约束)
% At_c 是 m×N² 行格式，SeDuMi 接受 N²×m 格式，故传入 At_c'
fprintf('\n=== 情形 2: SeDuMi (复数 SDP) ===\n');
K_sedumi          = struct();
K_sedumi.s        = N_c;
K_sedumi.scomplex = 1;
pars_sedumi.maxiter = 300;
pars_sedumi.fid     = 0;
obj_sed = NaN; t_sed = NaN;
maxc_c = max(abs(c_c));
tic;
try
    [x_sed, y_sed, info_sed] = sedumi(At_c', b_c, c_c/maxc_c, K_sedumi, pars_sedumi);
    t_sed = toc;
    if info_sed.numerr == 0
        fprintf('SeDuMi 成功，用时 %.2fs (迭代 %d 次)\n', t_sed, info_sed.iter);
    else
        fprintf('SeDuMi 数值问题 (numerr = %d)，用时 %.2fs\n', info_sed.numerr, t_sed);
    end
    obj_sed  = real(c_c' * x_sed);
    pinf_sed = norm(At_c * x_sed - b_c) / max(1, norm(b_c));
    S_sed    = reshape(c_c - maxc_c * At_c' * y_sed, N_c, N_c);
    S_sed    = (S_sed + S_sed') / 2;
    dS_sed   = real(eig(full(S_sed)));
    dinf_sed = max(0, -min(dS_sed)) / (1 + max(dS_sed));
    gap_sed  = abs(c_c'*x_sed - b_c'*y_sed*maxc_c) / (1 + abs(c_c'*x_sed) + abs(b_c'*y_sed*maxc_c));
    fprintf('  f = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', ...
            obj_sed, pinf_sed, dinf_sed, gap_sed);
catch ME
    fprintf('SeDuMi 出错: %s\n', ME.message);
end

%}
%% 5. convertCtoR 转换
fprintf('\n正在执行 CSDP → RSDP 转换...\n');
% 含对角约束（供 MOSEK 使用）
[At_r, b_r, c_r, K_r]             = convertCtoR(At_c,   b_c,   c_c,   K_c);
% 无对角约束（供 ManiSDP_complexunitdiag，由流形处理）
[At_r_nd, b_r_nd, c_r_nd, K_r_nd] = convertCtoR(At_cnd, b_cnd, c_cnd, K_cnd);
n_r = K_r_nd.s;

%{
%% 6. MOSEK (求解实 RSDP，含完整约束)
% At_r 是 (2N)²×m 列格式，convert_sedumi2mosek 直接接受该格式
fprintf('\n=== 情形 3: MOSEK (实 RSDP) ===\n');
max_time = 10000;
obj_mos = NaN; t_mosek = NaN;
maxc_r = max(abs(c_r));
try
    prob  = convert_sedumi2mosek(At_r, b_r, c_r/maxc_r, K_r);
    param.MSK_DPAR_OPTIMIZER_MAX_TIME = max_time - 10;
    param.MSK_IPAR_LOG = 0;
    tic;
    [rcode, res] = mosekopt('minimize echo(0)', prob, param);
    t_mosek = toc;
    if rcode == 0
        K_mosek.s = K_r.s;
        [X_mosek, y_mosek, S_mosek, mobj] = recover_mosek_sol_blk(res, K_mosek);
        if ~isempty(mobj)
            x_mos    = X_mosek{1}(:);
            obj_mos  = mobj(1) * maxc_r;
            pinf_mos = norm(At_r'*x_mos - b_r) / max(1, norm(b_r));
            by_mos   = b_r'*y_mosek * maxc_r;
            gap_mos  = abs(obj_mos - by_mos) / (1 + abs(obj_mos) + abs(by_mos));
            dS_mos   = eig(S_mosek{1});
            dinf_mos = max(0, -min(dS_mos)) / (1 + max(dS_mos));
            fprintf('MOSEK 成功，用时 %.2fs\n', t_mosek);
            fprintf('  f = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', ...
                    obj_mos, pinf_mos, dinf_mos, gap_mos);
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

%% 7. ManiSDP_complexunitdiag (求解实 RSDP)
fprintf('\n=== 情形 4: ManiSDP_complexunitdiag (实 RSDP) ===\n');
options_r.tol         = 1e-8;
options_r.AL_maxiter  = 1000;
options_r.p0          = 2;
options_r.sigma0      = 1e-3;
options_r.sigma_min   = 1e-2;
options_r.sigma_max   = 1e7;
options_r.gama        = 2;
options_r.tau1        = 1;
options_r.tau2        = 1;
options_r.theta       = 1e-3;
options_r.delta       = 8;
options_r.alpha       = 0.1;
options_r.TR_maxinner = 25;
options_r.TR_maxiter  = 4;
options_r.line_search = 0;

[X_r_mani, fval_r, data_r] = ManiSDP_complexunitdiag(At_r_nd, b_r_nd, c_r_nd, K_r_nd, options_r);
fprintf('ManiSDP_complexunitdiag: f = %.8f, eta = %.1e, t = %.2fs\n', ...
        fval_r, max([data_r.gap, data_r.pinf, data_r.dinf]), data_r.time);

 %{
%% 8. 结果对比
fprintf('\n=== 结果对比 ===\n');
fprintf('ManiCSDP_unitdiag:        f = %.8f, t = %.4fs\n', fval_c, data_c.time);
if ~isnan(obj_sed),  fprintf('SeDuMi (复数 SDP):        f = %.8f, t = %.4fs\n', obj_sed, t_sed); end
if ~isnan(obj_mos),  fprintf('MOSEK  (实 RSDP):         f = %.8f, t = %.4fs\n', obj_mos, t_mosek); end
fprintf('ManiSDP_complexunitdiag: f = %.8f, t = %.4fs\n', fval_r, data_r.time);
fprintf('[bqp_complex] script finished.\n');
drawnow;
%}
