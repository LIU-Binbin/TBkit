function C = gridvertcat(A, B)
    a = size(A, 2);
    b = size(B, 2);
    % 扩展A的列（整体重复b次），扩展B的列（每列重复a次）
    A_exp = repmat(A, 1, b);
    B_exp = repelem(B, 1, a);
    % 垂直拼接结果
    C = [A_exp; B_exp];
end