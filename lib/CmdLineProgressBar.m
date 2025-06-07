classdef CmdLineProgressBar < handle
    % CMDLINEPROGRESSBAR - 命令行进度条通知类
    % 用于在MATLAB命令窗口中显示动态更新的进度条
    %
    % 使用示例：
    %   pb = CmdLineProgressBar('Processing: '); 
    %   for k = 1 : 100
    %       pb.print(k, 100);  % 更新进度 (k/100)
    %       pause(0.02);       % 模拟耗时操作
    %   end
    %
    %   % 多进度条示例：
    %   pb = CmdLineProgressBar('Multi-task: ');
    %   for k = 1:10
    %       pb.print([k, 10-k], {['/10 '], ['/10 ']}); 
    %       pause(0.1);
    %   end
    %
    % 作者: Itamar Katz, itakatz@gmail.com
    % 修改: parkman
    % 优化: DeepSeek (修复析构函数递归问题，增强注释和健壮性)
    
    properties (Access = private)
        last_msg_len = 0;  % 存储上一次打印消息的长度（用于计算退格数量）
        is_completed = false; % 标记进度是否已完成（防止重复清理）
    end
    
    methods
        %% 构造函数 - 初始化进度条
        function obj = CmdLineProgressBar(msg)
            % 输入: 
            %   msg - 进度条前缀消息（字符串）
            fprintf('%s', msg);  % 打印初始消息（不换行）
        end
        
        %% 核心方法 - 更新进度显示
        function print(obj, n, tot, msg_tail)
            % 更新命令行进度显示
            %
            % 输入:
            %   n        - 当前进度值（标量或向量）
            %   tot      - 总进度值（标量或元胞数组）
            %   msg_tail - (可选) 附加在进度后的自定义消息
            %
            % 说明:
            %   - 单进度模式: n和tot为标量 (e.g., 5, 100)
            %   - 多进度模式: n为向量, tot为等长元胞数组 
            %        (e.g., n=[3,7], tot={'/10 ','/20 '})
            
            % 处理可选参数
            if nargin < 4
                msg_tail = '';  % 默认无附加消息
            end
            
            % 多进度条模式处理
            if numel(n) > 1 && numel(tot) == numel(n)
                obj.clearLastLine();  % 清除上一行显示
                
                % 构建多进度信息字符串 (e.g., "3/10 7/20")
                info_str = '';
                for i = 1:numel(n)
                    % 处理数值或字符串类型的进度单位
                    if iscell(tot)
                        unit_str = tot{i};
                    else
                        unit_str = num2str(tot(i));
                    end
                    info_str = [info_str, num2str(n(i)), unit_str];
                end
                info_str = [info_str, msg_tail];  % 附加尾部消息
                
                fprintf('%s', info_str);  % 打印新进度
                obj.last_msg_len = length(info_str);  % 更新长度记录
                
            % 单进度条模式处理
            else
                obj.clearLastLine();  % 清除上一行显示
                
                % 构建进度信息字符串
                if isscalar(n) && isscalar(tot)
                    info_str = sprintf('%d/%d', n, tot);  % 格式化输出
                else
                    info_str = [num2str(n), '/', num2str(tot)]; % 直接拼接
                end
                info_str = [info_str, msg_tail];  % 附加尾部消息
                
                fprintf('%s', info_str);
                obj.last_msg_len = length(info_str);  % 更新长度记录
                
                % 进度完成时标记并换行
                if n == tot
                    fprintf('\n');
                    obj.is_completed = true;  % 设置完成标志
                end
            end
        end
        
        %% 析构函数 - 清理进度显示
        function delete(obj)
            % 清理进度条显示（避免残留字符）
            if ~obj.is_completed && obj.last_msg_len > 0
                obj.clearLastLine();  % 清除进度显示
                fprintf('\n');        % 换行确保后续输出正常
            end
        end
    end
    
    methods (Access = private)
        %% 辅助方法 - 清除上一行显示
        function clearLastLine(obj)
            % 使用退格符清除上次打印的内容
            if obj.last_msg_len > 0
                fprintf('%s', char(8 * ones(1, obj.last_msg_len))); % ASCII 8 = 退格符
            end
        end
    end
end