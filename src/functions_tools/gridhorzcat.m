function C = gridhorzcat(A, B)
    m = size(A, 1);
    n = size(B, 1);
    % 创建行索引：A整体重复n次，B每行重复m次
    idxA = repmat((1:m)', n, 1);
    idxB = repelem((1:n)', m);
    % 直接索引组合结果
    C = [A(idxA, :), B(idxB, :)];
end