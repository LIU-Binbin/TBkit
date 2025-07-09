function [HR_list, eigval_list] = subsystem_construct(HR, Oper)
% SUBSYSTEM_CONSTRUCT 通过对称算符的本征空间构造多个子系统
%
%   [HR_list, eigval_list] = subsystem_construct(HR, Oper)
%
%   输入:
%       HR     - HR类对象（tight-binding 哈密顿量）
%       Oper   - Oper类对象，对称算符，含幺正矩阵 Oper.U
%
%   输出:
%       HR_list     - 一个 cell 数组，每个元素是一个子系统 HR 类
%       eigval_list - 一个向量，记录对应每个子系统的本征值（复数）

% 计算对称算符的本征分解
[V, D] = eig(Oper.U);
eigvals = diag(D);

% 数值上分类不同的本征值（去重，保留顺序）
tol = 1e-3;  % 容差
eigval_list = [];
group_idx = [];

for i = 1:length(eigvals)
    matched = false;
    for j = 1:length(eigval_list)
        if abs(eigvals(i) - eigval_list(j)) < tol
            group_idx(i) = j;
            matched = true;
            break;
        end
    end
    if ~matched
        eigval_list(end+1) = eigvals(i); %#ok<AGROW>
        group_idx(i) = length(eigval_list);
    end
end

% 显示本征值列表
fprintf('共识别出 %d 个子空间，对应本征值为：\n', length(eigval_list));
for k = 1:length(eigval_list)
    fprintf('  子系统 %d: 本征值 ≈ %s\n', k, num2str(eigval_list(k)));
end

% 为每个子系统构建 HR
HR_list = cell(1, length(eigval_list));

for k = 1:length(eigval_list)
    idx = (group_idx == k);
    Vk = V(:, idx);
    
    HRk = HR;
    HRk.HnumL = zeros(size(Vk,2), size(Vk,2), HR.NRPTS);
    for t = 1:HR.NRPTS
        HRk.HnumL(:,:,t) = Vk' * HR.HnumL(:,:,t) * Vk;
    end
    HRk.orbL = Vk' * HR.orbL;
    HRk.quantumL = Vk' * HR.quantumL;
    HRk.orb_symL = Vk' * HR.orb_symL;
    HR_list{k} = HRk;
end
end
