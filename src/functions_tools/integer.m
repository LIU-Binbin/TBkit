function A = integer(A)
    % 将矩阵A中的每个元素四舍五入到最近的整数
    % Round each element of matrix A to the nearest integer
    
    % 获取矩阵的大小
    % Get the size of the matrix
    [rows, cols] = size(A);  % 使用内置函数获取矩阵维度
    
    % 遍历矩阵的每个元素
    % Traverse each element in the matrix
    for i = 1:rows
        for j = 1:cols
            % 计算当前元素与最近整数的差值
            % Calculate the difference between current element and its nearest integer
            roundedValue = round(A(i,j));  % 先进行四舍五入
            
            % 判断是否需要替换
            % Determine whether to replace the value
            if abs(A(i,j) - roundedValue) < 1e-6
                A(i,j) = roundedValue;  % 替换为整数
            end
        end
    end
    
    % 返回处理后的矩阵
    % Return the processed matrix
end