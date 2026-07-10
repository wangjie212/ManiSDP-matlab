function [At_real, b_real, c_real, K_real] = convertCtoR(At, b, c, K)
    % 输入: At : m × n^2 (每行是一个向量化的约束矩阵)
    %       b  : m × 1
    %       c  : n^2 × 1
    % 输出: At_real : (2n)^2 × m 稀疏矩阵 (每列对应一个实约束)
    %       b_real  : m × 1
    %       c_real  : (2n)^2 × 1
    %       K_real  : 锥结构，K_real.s = 2n
    %
    % 转换关系: trace([A_R,-A_I;A_I,A_R] · X_r) = 2·trace(A_c · X_c)
    % 故约束矩阵除以 2，b 保持不变，c 同样除以 2。
    % 对 Hermitian A_c 和 X_c，Im(trace(A_c·X_c))=0 恒成立，无需虚部约束。

    n = K.s;
    n2 = 2 * n;
    Nreal = n2 * n2;
    m = length(b);

    % --- 目标向量 c_real ---
    C_mat = reshape(c, n, n);
    C_mat = (C_mat + C_mat') / 2;
    C_R = real(C_mat);
    C_I = imag(C_mat);
    C_real_mat = [C_R, -C_I; C_I, C_R];
    c_real = C_real_mat(:) / 2;

    % --- 逐约束构造三元组 ---
    row_idx = [];
    col_idx = [];
    vals    = [];
    b_real  = zeros(m, 1);

    for i = 1:m
        Ai_vec = full(At(i, :));
        Ai_mat = reshape(Ai_vec, n, n);
        Ai_R = real(Ai_mat);
        Ai_I = imag(Ai_mat);

        % 实嵌入：trace([A_R,-A_I;A_I,A_R]/2 · X_r) = trace(A_c · X_c) = b(i)
        Ai_real_mat = [Ai_R, -Ai_I; Ai_I, Ai_R] / 2;
        vec_real = Ai_real_mat(:);
        idx = find(vec_real);
        if ~isempty(idx)
            row_idx = [row_idx; idx];
            col_idx = [col_idx; i * ones(length(idx), 1)];
            vals    = [vals; vec_real(idx)];
        end
        b_real(i) = real(b(i));
    end

    % --- 构建最终稀疏矩阵 (Nreal × m) ---
    At_real = sparse(row_idx, col_idx, vals, Nreal, m);

    % --- 输出锥 ---
    K_real = K;
    K_real.s = n2;

    % 确保实数
    At_real = real(At_real);
    b_real  = real(b_real);
    c_real  = real(c_real);
end
