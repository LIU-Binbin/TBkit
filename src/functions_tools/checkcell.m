function [outputArg1] = checkcell(input1, cellobj)
% CHECKCELL 检查输入值是否存在于元胞数组中
%
% 输入参数：
%   input1 - 需要检查的值，可以是任意类型
%   cellobj - 元胞数组，包含多个元素
%
% 输出参数：
%   outputArg1 - 逻辑值，如果 input1 存在于 cellobj 中，则为 true；否则为 false
%
% 示例：
%   cellobj = {1, 'a', 3.14, [1, 2]};
%   outputArg1 = checkcell('a', cellobj); % 返回 true
%   outputArg1 = checkcell(2, cellobj);   % 返回 false

% 初始化输出为逻辑数组，默认值为 false
outputArg1 = false;

% 遍历元胞数组中的每个元素
for i = 1:length(cellobj)
    % 检查当前元胞元素是否与输入值相等
    if isequal(input1, cellobj{i})
        % 如果相等，将输出设置为 true 并退出循环
        outputArg1 = true;
        break;
    end
end
end

% function [outputArg1] = checkcell(input1,cellobj)
% %UNTITLED2 此处提供此函数的摘要
% %   此处提供详细说明
% outputArg1 = ~logical(1:length(cellobj));
% for i = 1:length(cellobj)
%     if  isequal(input1,cellobj{i})
%         outputArg1 = true;
%     end
% end
% end