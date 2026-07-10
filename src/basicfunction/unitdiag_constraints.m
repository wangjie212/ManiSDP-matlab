function [At, b, c, K] = unitdiag_constraints(C)
% 生成单位对角约束 SDP 的 SeDuMi 格式数据
%   min  ⟨C, X⟩
%   s.t. X ⪰ 0,  diag(X) = 1
% 输入:
%   C - n×n 实对称或复 Hermitian 矩阵（目标矩阵）
% 输出:
%   At - (n^2)×n 稀疏矩阵，每列对应一个约束 X_{ii}=1
%   b  - n×1 全1向量
%   c  - n^2×1 向量，C 按列拉直（复数）
%   K  - 结构体，K.s = n

    n = size(C, 1);
    % 目标函数系数
    c = C(:);               % 列向量化
    % 约束矩阵 At
    At = sparse(n^2, n);
    for i = 1:n
        At((i-1)*n + i, i) = 1;
    end
    b = ones(n, 1);
    K.s = n;
    K.l = 0;
end
