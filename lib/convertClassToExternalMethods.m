function convertClassToExternalMethods(classDefFile)
    % 验证输入是否为 .m 文件
    [classDir, className, ext] = fileparts(classDefFile);
    if ~strcmpi(ext, '.m')
        error('必须提供以 .m 结尾的类定义文件');
    end

    % 构造类文件夹路径 (@ClassName)
    classFolder = fullfile(classDir);
    privateFolder = fullfile(classFolder, 'private');

    % 创建文件夹（如果不存在）
    if ~exist(classFolder, 'dir')
        mkdir(classFolder);
    end
    if ~exist(privateFolder, 'dir')
        mkdir(privateFolder);
    end

    % 读取并处理类定义文件
    try
        classText = fileread(classDefFile);
    catch
        error('无法读取文件: %s', classDefFile);
    end

    [newClassDef, methodsToExtract] = processClassText(classText, className);

    % 写入更新后的类定义
    writeNewClassDef(classDefFile, newClassDef);

    % 将方法写入独立文件
    extractMethodsToFiles(methodsToExtract, classFolder, privateFolder);

    disp(['转换完成：方法已外部化为单独文件，请验证类 ' className ' 的功能是否正常。']);
end


function [newClassDef, methodsToExtract] = processClassText(classText, className)
    newClassDef = {};
    methodsToExtract = struct('name',{}, 'code',{}, 'isStatic',{}, 'isPrivate',{});
    
    % 提取 classdef 行
    classdefLine = regexp(classText, '^\s*classdef\s+[^\n\r]*', 'match', 'once', 'lineanchors');
    if isempty(classdefLine)
        error('无效的类定义：未找到 classdef 声明');
    end
    newClassDef{end+1} = strtrim(classdefLine);
    
    % 去掉 classdef 行后面的内容
    contentWithoutClassdef = regexprep(classText, '^\s*classdef\s+[^\n\r]*[\n\r]?', '', 'once', 'lineanchors');
    
    % 去掉注释、拆行
    lines = splitLinesAndRemoveComments(contentWithoutClassdef);
    
    % 逐行解析块结构
    currentBlockType = '';
    currentBlockLines = {};
    blockDepth = 0;

    % SavedBlock
    currentBlockTypeL = {};
    currentBlockLinesL = {};
    for i = 1:length(lines)
        rawLine = lines{i};
        line = strtrim(rawLine);

        if isempty(currentBlockType)
            blockStartMatch = regexpi(line, '^\s*(methods|properties|events)', 'once', 'match');
            if ~isempty(blockStartMatch)
                currentBlockType = blockStartMatch;
                blockDepth = 1;
                currentBlockLines = {rawLine};
                continue;
            end
        elseif ~isempty(currentBlockType)
            % 更新块深度
            pattern_tmp ='^\s*\<(methods|switch|arguments|properties|events|if|parfor|for|while|try|function)\>';
            if regexp(line, pattern_tmp, 'once')
                blockDepth = blockDepth + 1;
            elseif regexp(line, '^\s*end\s*;?\s*$', 'once')
                blockDepth = blockDepth - 1;
            end
            currentBlockLines{end+1} = rawLine;

            if blockDepth == 0
                currentBlockTypeL{end+1} = (currentBlockType);
                currentBlockLinesL{end+1} = currentBlockLines;
                currentBlockType = '';
                currentBlockLines = {};
            end
        else
            %newClassDef{end+1} = rawLine;
        end
    end

    for i =1:length(currentBlockLinesL)  
        [blockDecl, methods] = processBlock(currentBlockTypeL{i}, currentBlockLinesL{i}, className);
        newClassDef{end+1} = blockDecl;
        methodsToExtract = [methodsToExtract, methods];
    end

end


function lines = splitLinesAndRemoveComments(text)
    rawLines = splitlines(text);
    lines = {};
    inBlockComment = false;

    for i = 1:length(rawLines)
        line = strtrim(rawLines{i});

        % 检查是否是块注释起始（可能和其他内容混合）
        if inBlockComment
            if contains(line, '%}')
                inBlockComment = false;
            end
            continue;
        elseif contains(line, '%{')
            if ~contains(line, '%}')
                inBlockComment = true;
            end
            continue;
        end
       % 检查 fprintf 语句：如果是 fprintf 语句中的格式符号，则不当作注释处理
        if contains(line, 'fprintf')
            lines{end+1} = line; % 保留整行，避免误处理
            continue;
        end
        % 去除行尾注释，但保留 '...' 前的部分
        percentIdx = strfind(line, '%');
        if ~isempty(percentIdx)
            commentStart = percentIdx(1);
            codePart = strtrim(line(1:commentStart-1));
        else
            codePart = line;
        end

        if ~isempty(codePart)
            lines{end+1} = codePart;
        end
    end
end



function [name, code] = parseMethod(methodLines, className)
    code = strjoin(methodLines, '\n');
    name = '';
    if isempty(methodLines)
        return;
    end

    firstLine = strtrim(methodLines{1});
    
    % 正则匹配方法名：
    % 支持 function [out] = methodName(...) 或 function methodName(...)
    tokens = regexp(firstLine, ...
        '^\s*function\s+(?:\[?.*?\]?\s*=\s*)?(\w+)\s*\(', ...
        'tokens', 'once');

    if ~isempty(tokens)
        name = tokens{1};
        if strcmp(name, className)
            name = '';  % 构造函数不导出
        end
    end
end


function [blockDecl, methods] = processBlock(blockTypeLine, blockLines, className)
    methods = struct('name',{}, 'code',{}, 'isStatic',{}, 'isPrivate',{});
    
    % 判断是否 methods 块
    if ~startsWith(strtrim(lower(blockTypeLine)), 'methods')
        blockDecl = strjoin(blockLines, '\n');
        return;
    end

    modifiers = parseMethodModifiers(blockTypeLine);
    [methodSections, declarations] = extractMethodSections(blockLines);

    % 构建 methods 声明
    if ~isempty(declarations)
        blockDecl = sprintf('%s',  strjoin(declarations, '\n'));
    else
        blockDecl = strtrim(blockTypeLine);
    end

    % 提取方法名和代码
    for i = 1:length(methodSections)
        [name, code] = parseMethod(methodSections{i}, className);
        if ~isempty(name)
            methods(end+1) = struct( ...
                'name', name, ...
                'code', code, ...
                'isStatic', modifiers.isStatic, ...
                'isPrivate', modifiers.isPrivate ...
            );
        end
    end
end



function modifiers = parseMethodModifiers(blockLines)
    modifiers.isStatic   = false;
    modifiers.isPrivate  = false;
    modifiers.access     = 'public';
    modifiers.isAbstract = false;
    modifiers.isSealed   = false;
    modifiers.isHidden   = false;
    
    % 拼接声明行并去掉续行符
    fullDecl = regexprep(join(blockLines, ' '), '\.\.\.\s*', '');

    % 检测常见修饰符
    if contains(fullDecl, 'Static', 'IgnoreCase', true)
        modifiers.isStatic = true;
    end
    if contains(fullDecl, 'Abstract', 'IgnoreCase', true)
        modifiers.isAbstract = true;
    end
    if contains(fullDecl, 'Sealed', 'IgnoreCase', true)
        modifiers.isSealed = true;
    end
    if contains(fullDecl, 'Hidden', 'IgnoreCase', true)
        modifiers.isHidden = true;
    end
    
    % Access 级别
    accessMatch = regexpi(fullDecl, 'Access\s*=\s*(\w+)', 'tokens', 'once');
    if ~isempty(accessMatch)
        modifiers.access = lower(accessMatch{1});
        modifiers.isPrivate = strcmpi(modifiers.access, 'private');
    end
end




function [methodSections, declarations] = extractMethodSections(blockLines)
    methodSections = {};
    declarations = {};
    inMethod = false;
    currentMethod = {};
    depth = 0;
    
    for i = 1:length(blockLines)
        line = strtrim(blockLines{i});
        rawLine = blockLines{i};
        
        if regexp(line, '^\s*\<(function)\>', 'once')
            if ~inMethod
                inMethod = true;
                depth = 1;
                currentMethod = {rawLine};
                declLine = regexprep(rawLine, '\s*$', ';');
                declLine = regexprep(rawLine, 'function', '');
                declarations{end+1} = declLine;
            else
                currentMethod{end+1} = rawLine;
                depth = depth + 1;
            end
        elseif inMethod
            currentMethod{end+1} = rawLine;
            if regexp(line, '^\s*end\s*;?\s*$', 'once')
                depth = depth - 1;
                if depth == 0
                    inMethod = false;
                    methodSections{end+1} = currentMethod;
                    currentMethod = {};
                end
            elseif regexp(line, '^\s*\<(switch|arguments|properties|events|if|parfor|for|while|try|function)\>', 'once')
                depth = depth + 1;
                % disp(depth);
            end
        else
            declarations{end+1} = rawLine;
        end
    end

    if inMethod
        warning('有方法可能缺少 `end` 结尾，请检查方法定义：\n%s', currentMethod{1});
    end
end


function extractMethodsToFiles(methods, classFolder, privateFolder)
    for k = 1:length(methods)
        method = methods(k);
        destFolder = classFolder;
        if method.isPrivate
            destFolder = privateFolder;
        end

        if ~exist(destFolder, 'dir')
            mkdir(destFolder);
        end
        
        filePath = fullfile(destFolder, [method.name '.m']);
        fid = fopen(filePath, 'w');
        if fid == -1
            error('无法写入方法文件：%s', filePath);
        end
        fprintf(fid, '%s\n', method.code);
        fclose(fid);

        fprintf('写入方法: %s\n', filePath);
    end
end


function writeNewClassDef(filePath, classDef)
    % 检查文件夹是否存在，如果不存在则创建
    [folderPath, ~, ~] = fileparts(filePath);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    % 尝试打开文件
    fid = fopen(filePath, 'w');
    if fid == -1
        error('无法写入文件: %s', filePath);
    end
    
    % 写入类定义内容
    fprintf(fid, '%s\n', classDef{:});
    % 写入最后的end
    fprintf(fid, '%s\n', 'end');
    % 关闭文件
    fclose(fid);
    
    fprintf('类定义已成功写入: %s\n', filePath);
end
