
clear; 
clc;

rng(1);

d = 15;
N = (2*d+1)*(d+1);

% 生成随机四次型目标矩阵 Q
Q_quartic = randn(N) + 1i*randn(N);
Q_quartic = (Q_quartic + Q_quartic') / 2;

%% 调用 qsmom_complex 生成二阶松弛参数
% qsmom_complex 返回 A_sdp: m × n² (行格式), b_sdp, c_sdp, K
[A_sdp, b_sdp, c_sdp, K] = qsmom_complex(Q_quartic);

%% ManiCSDP
fprintf('\n--- ManiCSDP 正在求解 ---\n');
options.sigma0 = 1;
options.AL_maxiter = 500;
options.p0=1;
options.tol = 1e-8;
options.tau1=1e-2;
options.tau2=1e-2;
options.delta=15;
options.TR_maxinner=20;
options.TR_maxiter=6;
options.line_search = 1;
tic;
[X_mani, obj_mani, data_mani] = ManiCSDP(A_sdp', b_sdp, c_sdp, K, options);
emani = max([data_mani.gap, data_mani.pinf, data_mani.dinf]);
tmani = toc;
fprintf('ManiCSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', obj_mani, emani, tmani);


%{
%% SeDuMi
fprintf('调用 SeDuMi (复数 SDP, 最大迭代 200 次)\n');
K_sedumi = struct();
K_sedumi.s = K.s;       % n×n SDP 块大小
K_sedumi.scomplex = 1;  % 第 1 个块为复 Hermitian
pars.maxiter = 200;
pars.fid = 0;
tic;
try
    [x_sed, y_sed, info_sed] = sedumi(A_sdp', b_sdp, c_sdp, K_sedumi, pars);
    t_sedumi = toc;
    if info_sed.numerr == 0
        fprintf('SeDuMi 成功完成，用时 %.2f 秒 (迭代 %d 次)\n', t_sedumi, info_sed.iter);
    else
        fprintf('SeDuMi 数值问题 (numerr = %d)，用时 %.2f 秒，当前残差：\n', info_sed.numerr, t_sedumi);
    end
    obj_sed = real(c_sdp' * x_sed);
    pinf = norm(A_sdp*x_sed - b_sdp) / max(1, norm(b_sdp));
    S_mat = reshape(c_sdp - A_sdp' * y_sed, K.s, K.s);
    dS = eig((S_mat + S_mat') / 2);
    dinf = max(0, -min(real(dS))) / (1 + max(real(dS)));
    gap = abs(c_sdp'*x_sed - b_sdp'*y_sed) / (1 + abs(c_sdp'*x_sed) + abs(b_sdp'*y_sed));
    fprintf('  最优值 = %.8f, pinf = %.2e, dinf = %.2e, gap = %.2e\n', obj_sed, pinf, dinf, gap);
catch ME
    fprintf('SeDuMi 出错: %s\n', ME.message);
end
%}


%% convertCtoR + ManiSDP + MOSEK
% convertCtoR 接受 m × n² 行格式，返回 At_r: (2n)² × actual_m 列格式
[At_r, b_r, c_r, K_r] = convertCtoR(A_sdp, b_sdp, c_sdp, K);
n_r = K_r.s;


fprintf('\n--- ManiSDP 正在求解 ---\n');
options_r = struct();
options_r.p0 = 1;
options_r.sigma0 = 1;
options_r.sigma_min = 1e-1;
options_r.theta=1e-2;
options_r.delta=15;
options_r.tau1=1e-2;
options_r.tau2=1e-2;
options_r.TR_maxinner=20;
options_r.TR_maxiter=6;
options_r.line_search = 1;
tic;
[x_r, obj_r, data_r] = ManiSDP(At_r, b_r, c_r, K_r, options_r);
emani_r = max([data_r.gap, data_r.pinf, data_r.dinf]);
t_manisdp = toc;
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', obj_r, emani_r, t_manisdp);

%{
fprintf('调用 MOSEK (实数 SDP)\n');
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
            x_mos = X_mosek{1}(:);
            pinf = norm(At_r'*x_mos - b_r) / max(1, norm(b_r));
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

