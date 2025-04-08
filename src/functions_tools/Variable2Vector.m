function Vector = Variable2Vector(Variable)
    % Variable2Vector converts a variable name in the format 'prefix_mX_mY_mZ' 
    % into a numerical vector [X, Y, Z].
    %
    % Input:
    %   Variable - A string representing the variable name in the format 
    %              'prefix_mX_mY_mZ'.
    %
    % Output:
    %   Vector   - A 1x3 numerical vector [X, Y, Z].

    % Split the input string at underscores ('_') into a cell array of strings.
    % The first element is the prefix, and the subsequent elements are 'mX', 'mY', and 'mZ'.
    StrL = strsplit(Variable, '_');
    
    % Extract the elements corresponding to 'mX', 'mY', and 'mZ'.
    % These are the 2nd, 3rd, and 4th elements of the split string array.
    StrL = StrL(2:4);
    
    % Replace the character 'm' with a minus sign ('-') in each of the extracted strings.
    % This is done to handle negative values represented by 'mX', 'mY', and 'mZ'.
    StrL = strrep(StrL, 'm', '-');
    
    % Convert the modified strings to numerical values and store them in the output vector.
    Vector = str2double(StrL);
end
