function [A, b, c, K] = bqpmom_complex(Q, e)
    % 输入: Q (n x n Hermitian矩阵), e (n x 1 复向量)
    % 输出: A (m x N^2), b (m x 1), c (N^2 x 1), K.s = N
    
    n = length(e);
    
    % --- 第一步：构造基底并向量化 (2n维列表) ---
    % 每个基底向量表示为 [z_1...z_n; conj(z_1)...conj(z_n)]
    basis_vecs = generate_basis_vectors(n);
    N = length(basis_vecs);
    K.s = N; 
    
    % --- 第二步：从基底指标建立矩阵 X 各个位置的指标映射 ---
    % class_map: Key 是 n 维净指标向量的字符串, Value 是对应的矩阵线性索引数组
    class_map = containers.Map();
    
    for j = 1:N
        for i = 1:N
            % X(i,j) 对应的单项式向量 = basis{i} + conj(basis{j})
            % conj 操作对应向量化表示中：前 n 维与后 n 维互换
            vec_i = basis_vecs{i};
            vec_j_conj = flip_half_vector(basis_vecs{j}, n);
            combined_2n = vec_i + vec_j_conj;
            
            % --- 第三步：BQP 简化 (前一半减后一半) ---
            % 得到 n 维新列表（净指标），作为等价类的标签
            net_idx = combined_2n(1:n) - combined_2n(n+1:end);
            label = vec2str(net_idx);
            
            flat_idx = (j-1)*N + i;
            if isKey(class_map, label)
                class_map(label) = [class_map(label), flat_idx];
            else
                class_map(label) = flat_idx;
            end
        end
    end
    
    % --- 第四步：构造约束 A(X) = b ---
    At_list = {}; % 使用 cell 数组暂存列向量，最后统一拼接转置
    b = [];
    
    % 1. 显式处理常数项类 (净指标为 0 的类)，设定基准 X(1,1) = 1
    zero_label = vec2str(zeros(n, 1));
    if isKey(class_map, zero_label)
        zero_idxs = class_map(zero_label);
        
        % 约束 Re(X(1,1)) = 1
        [Ar, Ai] = get_single_pos_matrix(1, N);
        At_list{end+1} = sparse(Ar(:)); b = [b; 1];
        % 约束 Im(X(1,1)) = 0 (保证对角线实数性)
        At_list{end+1} = sparse(Ai(:)); b = [b; 0];
        
        % 让该类中其他位置等于 X(1,1)
        base_idx = zero_idxs(1);
        for m = 2:length(zero_idxs)
            [Ar, Ai] = get_diff_matrices(base_idx, zero_idxs(m), N);
            At_list{end+1} = sparse(Ar(:)); b = [b; 0];
            At_list{end+1} = sparse(Ai(:)); b = [b; 0];
        end
    end
    
    % 2. 处理其他“非零”类 (采用“砍半”逻辑，避免共轭冗余)
    all_labels = keys(class_map);
    processed_labels = containers.Map(); 
    processed_labels(zero_label) = true; % 跳过常数项
    
    for k = 1:length(all_labels)
        curr_label = all_labels{k};
        if isKey(processed_labels, curr_label), continue; end
        
        % 获取当前类及其共轭类的标签 (净指标取相反数即为共轭)
        curr_v = str2vec(curr_label);
        conj_label = vec2str(-curr_v); 
        
        % 标记当前类及其共轭类已处理
        processed_labels(curr_label) = true;
        processed_labels(conj_label) = true;
        
        % 类内元素相等约束
        idxs = class_map(curr_label);
        base_idx = idxs(1);
        for m = 2:length(idxs)
            [Ar, Ai] = get_diff_matrices(base_idx, idxs(m), N);
            At_list{end+1} = sparse(Ar(:)); b = [b; 0];
            At_list{end+1} = sparse(Ai(:)); b = [b; 0];
        end
    end
    
    % 统一转置，生成 ManiCSDP 需要的 m x N^2 矩阵
    A = cell2mat(At_list)';
    
    % --- 第五步：构造目标矩阵 C ---
    C_tmp = zeros(N, N);
    
    % 1. 映射一次项 e'z (净指标为 e_i 位置为 1)
    for i = 1:n
        v_e = zeros(n, 1); v_e(i) = 1;
        l_e = vec2str(v_e);
        if isKey(class_map, l_e)
            temp_idxs = class_map(l_e); % 修正：拆分链式调用
            target_idx = temp_idxs(1); 
            C_tmp(target_idx) = C_tmp(target_idx) + e(i);
        end
    end
    
    % 2. 映射二次项 z'Qz (净指标为 i 位 1, j 位 -1)
    for i = 1:n
        for j = 1:n
            v_Q = zeros(n, 1); v_Q(i) = 1; v_Q(j) = v_Q(j) - 1;
            l_Q = vec2str(v_Q);
            if isKey(class_map, l_Q)
                temp_idxs = class_map(l_Q); % 修正：拆分链式调用
                target_idx = temp_idxs(1);
                C_tmp(target_idx) = C_tmp(target_idx) + Q(i,j);
            end
        end
    end
    
    % 埃尔米特化并向量化
    C_mat = (C_tmp + C_tmp') / 2;
    c = C_mat(:);
end

%% --- 内部工具函数定义 ---

function vecs = generate_basis_vectors(n)
    % 生成二阶基底向量 [1; z; conj(z); z*z ...]
    vecs = {zeros(2*n, 1)}; % 常数 1
    for i = 1:n % z_i
        v = zeros(2*n, 1); v(i) = 1; vecs{end+1} = v;
    end
    %for i = 1:n % conj(z_i)
    %    v = zeros(2*n, 1); v(n+i) = 1; vecs{end+1} = v;
    %end
    % 此处根据需要可继续添加 z_i*z_j 等更高阶组合
    % zz
    for i = 1:n, for j = i:n, v = zeros(2*n, 1); v(i)=v(i)+1; v(j)=v(j)+1; vecs{end+1}=v; end; end
    % zconj
    %for i = 1:n, for j = 1:n, v = zeros(2*n, 1); v(i)=1; v(n+j)=1; vecs{end+1}=v; end; end
    % conjconj
    %for i = 1:n, for j = i:n, v = zeros(2*n, 1); v(n+i)=v(n+i)+1; v(n+j)=v(n+j)+1; vecs{end+1}=v; end; end
end

function v_out = flip_half_vector(v_in, n)
    % 对应 conj 操作：[z项; conj项] -> [conj项; z项]
    v_out = [v_in(n+1:end); v_in(1:n)];
end

function [Ar, Ai] = get_diff_matrices(idx1, idx2, N)
    % 构造 Re(X1 - X2) = 0 和 Im(X1 - X2) = 0 的对称矩阵
    A_tmp_r = zeros(N, N);
    A_tmp_r(idx1) = 1; A_tmp_r(idx2) = -1;
    Ar = (A_tmp_r + A_tmp_r') / 2;
    
    A_tmp_i = zeros(N, N);
    A_tmp_i(idx1) = 1i; A_tmp_i(idx2) = -1i;
    Ai = (A_tmp_i + A_tmp_i') / 2;
end

function [Ar, Ai] = get_single_pos_matrix(idx, N)
    % 构造提取单个位置实部/虚部的对称矩阵
    A = zeros(N, N); A(idx) = 1;
    Ar = (A + A') / 2;
    Ai = (A - A') / (2i);
end

function s = vec2str(v)
    s = sprintf('%d,', v);
end

function v = str2vec(s)
    v = str2num(strrep(s, ',', ' '))';
end