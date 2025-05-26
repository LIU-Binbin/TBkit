function latexCode = table2latex(dataTable)
    % 检查行名是否存在
    hasRowNames = ~isempty(dataTable.Properties.RowNames);
    [numRows, numVars] = size(dataTable);
    
    % 处理行名
    if hasRowNames
        rowNames = dataTable.Properties.RowNames;
        rowNameCells = cellfun(@escapeLatex, rowNames, 'UniformOutput', false);
    end
    
    % 处理列名
    varNames = dataTable.Properties.VariableNames;
    %varNameCells = cellfun(@escapeLatex, varNames, 'UniformOutput', false);
    
    % 转换数据为字符串并转义特殊字符
    dataCells = cell(numRows, numVars);
    for i = 1:numRows
        for j = 1:numVars
            value = dataTable{i, j};
            dataCells{i,j} = escapeLatex(value);
        end
    end
    
    % 构建列格式
    if hasRowNames
        colFormat = ['l', repmat('c', 1, numVars)];
    else
        colFormat = repmat('c', 1, numVars);
    end
    
    % 生成LaTeX代码
    latexLines = {};
    latexLines{end+1} = '\begin{tabular}';
    latexLines{end+1} = ['{', colFormat, '}'];
    latexLines{end+1} = '\hline';
    
    % 表头
    if hasRowNames
        header = [' ', varNames];
    else
        header = varNames;
    end
    headerStr = strjoin(header, ' & ');
    latexLines{end+1} = [headerStr, ' \\'];
    latexLines{end+1} = '\hline';
    
    % 数据行
    for i = 1:numRows
        if hasRowNames
            row = [rowNameCells(i), dataCells(i, :)];
        else
            row = dataCells(i, :);
        end
        rowStr = strjoin(row, ' & ');
        latexLines{end+1} = [rowStr, ' \\'];
    end
    
    latexLines{end+1} = '\hline';
    latexLines{end+1} = '\end{tabular}';
    
    % 合并为字符串
    latexCode = strjoin(latexLines, '\n');
end

% 转义LaTeX特殊字符
function str = escapeLatex(s)
if isstring(s)
    str = (string(s));
elseif isnumeric(s)
    str = (string(s));
elseif isa(s,'sym')
    str = latex(s);
    str = char("$"+str+"$");
else
    str = string(s);
end
    str = strrep(str,' ','');
    str = strrep(str,',','');
    str = strrep(str,'\\chi','\chi');
    %
    % replacements = {'\', '\\'; '&', '\&'; '_', '\_'; '%', '\%';
    %                 '$', '\$'; '#', '\#'; '{', '\{'; '}', '\}';
    %                 '^', '\^{}'; '~', '\~{}'};
    % for i = 1:size(replacements, 1)
    %     str = strrep(str, replacements{i,1}, replacements{i,2});
    % end
    %
end