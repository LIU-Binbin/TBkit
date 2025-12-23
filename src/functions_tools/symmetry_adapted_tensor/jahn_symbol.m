function Tensor = jahn_symbol(jahn_symbol_Str)
    % jahn_symbol_Char = char(jahn_symbol_Str);
    jahn_symbol_Str = string(jahn_symbol_Str);
    if strcontain(jahn_symbol_Str,"e")
        Tensor.is_pesudo_tensor = true;
    else
        Tensor.is_pesudo_tensor = false;
    end
    jahn_symbol_Str = strrep(jahn_symbol_Str,"e","");
    
    if strcontain(jahn_symbol_Str,"a")
        Tensor.is_Time_reversal_odd = true;
    else
        Tensor.is_Time_reversal_odd = false;
    end
    jahn_symbol_Str = strrep(jahn_symbol_Str,"a","");
    
    
    % Decompose
    % [[V2]3] into [V2][V2][V2]
    % omit temp
    
    %  [V2]{V2}[V2]' into [VV]{VV}[VV]'
    outputStr = regexprep(jahn_symbol_Str, '([A-Za-z])(\d+)', '${repmat($1, 1, str2double($2))}');
    Tensor.tensor_rank = sum(isletter(outputStr));
    Tensor.jahn  = outputStr;
    %
    % fprintf('%s : %d\n',outputStr,Tensor.tensor_rank);
    Tensor = CheckPair(Tensor);
    Tensor.n_ele = 3^Tensor.tensor_rank ;
    %
    Tensor = init(Tensor);
    end
    
    function Tensor = CheckPair(Tensor)
    % pair
    % Tensor.index = 1:Tensor.tensor_rank;
    Tensor.SymmetricPair = [];
    Tensor.AsymmetricPair = [];
    Tensor.TsymmetricPair = [];
    Tensor.TasymmetricPair = [];
    %
    CharJahn = char(Tensor.jahn);
    Pos = 0;
    StackFlag = 0;
    for ichar = 1:numel(CharJahn)
        CharTmp = CharJahn(ichar);
        if ichar == numel(CharJahn)
            ichar_plus1 = ichar;
        else
            ichar_plus1 = ichar+1;
        end
        switch CharTmp
            case {'[','{'}
                StackFlag = Pos + 1;
            case {']'}
    
                if strcmp(CharJahn(ichar_plus1),"*")
                    Tensor.TsymmetricPair = [Tensor.TsymmetricPair;StackFlag,Pos];
                else
                    Tensor.SymmetricPair = [Tensor.SymmetricPair;StackFlag,Pos];
                end
                StackFlag = 0;
            case {'}'}
                if strcmp(CharJahn(ichar_plus1),"*")
                    Tensor.TasymmetricPair = [Tensor.TasymmetricPair ;StackFlag,Pos];
                else
                    Tensor.AsymmetricPair = [Tensor.AsymmetricPair;StackFlag,Pos];
                end
                StackFlag = 0;
            case {'V'}
                Pos = Pos + 1;
            otherwise
                if isletter(CharTmp)
                    Pos = Pos + 1;
                else
    
                end
        end
    end
end

% %  bad
% inputStr = jahn_symbol_Str;
% % Step 0:
% outputStr = regexprep(inputStr, '\]', ')');
% outputStr = regexprep(outputStr, '\}', ')');
% outputStr = regexprep(outputStr, '\[', '(');
% outputStr = regexprep(outputStr, '\{', '(');
% % Step 1: 字母之间插入 *
% outputStr = regexprep(outputStr, '([A-Za-z])([A-Za-z])', '$1*$2');
%
% % Step 2: 字母后跟数字插入 ^
% outputStr = regexprep(outputStr, '([A-Za-z])(\d)', '$1^$2');
%
% % Step 3: 数字后跟字母或正括号插入 *
% outputStr = regexprep(outputStr, '(\d)([A-Za-z(])', '$1*$2');
%
% % Step 4: 反括号与正括号之间插入 *
% outputStr = regexprep(outputStr, '\)\(' , ')*(');
%
% %
% outputSym = simplify(str2sym(outputStr));
%  child  = children(V^6);
% Tensor.tensor_rank = double(child{2});
