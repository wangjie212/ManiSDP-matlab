function [A, b, c, K] = qsmom_complex(Q_quartic)
    % 1. 自动获取变量个数 n 并构造基底
    N_total = size(Q_quartic, 1);
    n = round((-3 + sqrt(9 - 8*(1-N_total))) / 4); % 由 N = (2n+1)(n+1) 反推 n
    
    basis_vecs = generate_basis_vectors_qs(n); 
    N = length(basis_vecs); 
    K.s = N;
    
    % 2. 建立全量等价类映射 (借鉴 BQP 的 Net-Label 逻辑)
    class_map = containers.Map();
    for j = 1:N
        for i = 1:N
            % 计算 X(i,j) 对应的单项式标签: basis{i} + conj(basis{j})
            % conj 对应前后 n 位互换 (Block Shift)
            v_i = basis_vecs{i};
            v_j_conj = [basis_vecs{j}(n+1:2*n); basis_vecs{j}(1:n)];
            
            label_vec = v_i + v_j_conj;
            label_str = sprintf('%d,', label_vec);
            
            flat_idx = (j-1)*N + i;
            if isKey(class_map, label_str)
                class_map(label_str) = [class_map(label_str), flat_idx];
            else
                class_map(label_str) = flat_idx;
            end
        end
    end
    
    At_list = {}; b_vals = [];

    % 3. [完整基底专用] 结构性约束 (Hankel 结构) - 同标签元素相等
    all_keys = keys(class_map);
    processed_labels = containers.Map();

    for k = 1:length(all_keys)
        curr_l = all_keys{k};
        if isKey(processed_labels, curr_l), continue; end

        % 寻找共轭标签
        v = sscanf(curr_l, '%d,');
        conj_l = sprintf('%d,', [v(n+1:2*n); v(1:n)]);

        processed_labels(curr_l) = true;
        processed_labels(conj_l) = true;

        idxs = class_map(curr_l);
        if length(idxs) > 1
            base_idx = idxs(1);
            for m = 2:length(idxs)
                [Ar, Ai] = get_diff_matrices(base_idx, idxs(m), N);
                add_to_A(Ar, 0); add_to_A(Ai, 0);
            end
        end
    end
    
    % 4. 球面定位约束 w_k * (sum |z_i|^2 - 1) = 0
    % Lasserre d=2 中定位矩阵 M_1(g·y) 由 deg ≤ 1 的基底元素索引：
    %   简化基底（仅 z 单项式）：{1, z_1,...,z_n}，共 1+n 个 → upper_k = 1+n
    %   完整基底（含 conj(z_i)）：{1, z_i, conj(z_i)}，共 1+2n 个 → upper_k = 1+2*n
    % 注：本循环只覆盖 M_1 的第一列 (h, 1) 即 case 1/2；(z_j, z_k) 块需另写双循环。
    %upper_k = 1 + n;          % 简化基底
    upper_k = 1 + 2*n;       % 完整基底

    % [完整基底专用] 共轭对去重：conj(w_k) 的约束 = w_k 约束的共轭，Hermitian PSD 自动满足
    processed_sphere = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    for k = 1:upper_k
        v_wk = basis_vecs{k};

        % [完整基底专用] 跳过共轭冗余
        wk_label   = sprintf('%d,', v_wk);
        conj_label = sprintf('%d,', [v_wk(n+1:2*n); v_wk(1:n)]);
        if isKey(processed_sphere, conj_label) && ~strcmp(wk_label, conj_label)
            continue;
        end
        processed_sphere(wk_label) = true;

        M_re = sparse(N, N); M_im = sparse(N, N);

        % (A) 累加 w_k * |z_i|^2
        for i = 1:n
            v_target = v_wk;
            v_target(i)   = v_target(i)   + 1;
            v_target(n+i) = v_target(n+i) + 1;

            t_label = sprintf('%d,', v_target);
            if isKey(class_map, t_label)
                idx_list = class_map(t_label);
                [Ar, Ai] = get_single_pos_matrix(idx_list(1), N);
                M_re = M_re + Ar; M_im = M_im + Ai;
            end
        end

        if k == 1
            % 0 阶球面约束：sum_i |z_i|^2 = 1（直接 b=1，不耦合 X(1,1)；M_im=0 自动跳过）
            add_to_A(M_re, 1);
        else
            % (B) 减去 w_k * 1，得 sum_i w_k|z_i|^2 - w_k = 0
            t_label_wk = sprintf('%d,', v_wk);
            idx_wk = class_map(t_label_wk);
            [Ar_wk, Ai_wk] = get_single_pos_matrix(idx_wk(1), N);

            add_to_A(M_re - Ar_wk, 0);
            add_to_A(M_im - Ai_wk, 0);
        end
    end
    %{
    % 原始实现（含冗余共轭约束，备用恢复）
    for k = 1:N
        v_wk = basis_vecs{k};
        M_re = sparse(N, N); M_im = sparse(N, N);
        for i = 1:n
            v_target = v_wk;
            v_target(i) = v_target(i) + 1;
            v_target(n+i) = v_target(n+i) + 1;
            t_label = sprintf('%d,', v_target);
            if isKey(class_map, t_label)
                idx_list = class_map(t_label);
                [Ar, Ai] = get_single_pos_matrix(idx_list(1), N);
                M_re = M_re + Ar; M_im = M_im + Ai;
            end
        end
        t_label_wk = sprintf('%d,', v_wk);
        idx_wk = class_map(t_label_wk);
        [Ar_wk, Ai_wk] = get_single_pos_matrix(idx_wk(1), N);
        add_to_A(M_re - Ar_wk, 0);
        add_to_A(M_im - Ai_wk, 0);
    end
    %}

    % 4b. (z_j, z_k) 块定位约束 z_j·conj(z_k)·(sum|z_i|^2 - 1) = 0
    %     即 sum_i X(idx(z_j z_i), idx(z_k z_i)) = X(1+j, 1+k)，对所有 j ≤ k
    %     j = k：1 个实约束；j < k：实部 + 虚部各 1 个
    %     对应定位矩阵 M_1(g·y) 的 (j+1, k+1) 块，是 case 1/2 之外的剩余约束
    for j = 1:n
        for k = j:n
            M_re = sparse(N, N); M_im = sparse(N, N);

            % (A) 累加 sum_i E[z_j z_i · conj(z_k z_i)]
            for i = 1:n
                v_target = zeros(2*n, 1);
                v_target(j)   = v_target(j)   + 1;
                v_target(i)   = v_target(i)   + 1;
                v_target(n+k) = v_target(n+k) + 1;
                v_target(n+i) = v_target(n+i) + 1;

                t_label = sprintf('%d,', v_target);
                if isKey(class_map, t_label)
                    idx_list = class_map(t_label);
                    [Ar, Ai] = get_single_pos_matrix(idx_list(1), N);
                    M_re = M_re + Ar; M_im = M_im + Ai;
                end
            end

            % (B) 减去 E[z_j conj(z_k)]
            v_jk = zeros(2*n, 1);
            v_jk(j)   = 1;
            v_jk(n+k) = 1;
            t_label_jk = sprintf('%d,', v_jk);
            idx_jk = class_map(t_label_jk);
            [Ar_jk, Ai_jk] = get_single_pos_matrix(idx_jk(1), N);

            add_to_A(M_re - Ar_jk, 0);
            add_to_A(M_im - Ai_jk, 0);  % j = k 时差为零矩阵，add_to_A 自动跳过
        end
    end

    % 5. 显式设定 X(1,1) = 1 (常数项基准)
    [Ar11, Ai11] = get_single_pos_matrix(1, N);
    add_to_A(Ar11, 1); add_to_A(Ai11, 0);

    % 6. 输出结果
    A = cell2mat(At_list)';
    b = b_vals;
    c = Q_quartic(:);

    % ================== 嵌套辅助函数 ==================
    
    function add_to_A(M, val)
        if nnz(M) > 0
            At_list{end+1} = sparse(M(:));
            b_vals(end+1, 1) = val;
        end
    end

    function [Ar, Ai] = get_diff_matrices(idx1, idx2, dim)
        % 构造 Re(X1 - X2) = 0 和 Im(X1 - X2) = 0
        A_tmp_r = zeros(dim, dim);
        A_tmp_r(idx1) = 1; A_tmp_r(idx2) = -1;
        Ar = (A_tmp_r + A_tmp_r') / 2;
        
        A_tmp_i = zeros(dim, dim);
        A_tmp_i(idx1) = 1i; A_tmp_i(idx2) = -1i;
        Ai = (A_tmp_i + A_tmp_i') / 2;
    end

    function [Ar, Ai] = get_single_pos_matrix(idx, dim)
        % 提取单个位置的实部和虚部
        A = zeros(dim, dim); A(idx) = 1;
        Ar = (A + A') / 2;
        Ai = (A - A') / (2i);
    end
end

% ================== 外部子函数 ==================

function vecs = generate_basis_vectors_qs(n)
    % 严格按照 [1; z; conj(z); zz; zconj; conjconj] 顺序
    vecs = {zeros(2*n, 1)}; % 1
    for i = 1:n, v = zeros(2*n, 1); v(i)=1; vecs{end+1}=v; end % z
    for i = 1:n, v = zeros(2*n, 1); v(n+i)=1; vecs{end+1}=v; end % conj(z)
    % zz
    for i = 1:n, for j = i:n, v = zeros(2*n, 1); v(i)=v(i)+1; v(j)=v(j)+1; vecs{end+1}=v; end; end
    % zconj
    for i = 1:n, for j = 1:n, v = zeros(2*n, 1); v(i)=1; v(n+j)=1; vecs{end+1}=v; end; end
    % conjconj
    for i = 1:n, for j = i:n, v = zeros(2*n, 1); v(n+i)=v(n+i)+1; v(n+j)=v(n+j)+1; vecs{end+1}=v; end; end
end