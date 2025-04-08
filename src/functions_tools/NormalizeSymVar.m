function SymVarList = NormalizeSymVar(SymVarList)
    % NormalizeSymVar processes a list of symbolic variables, breaking them into
    % real and imaginary parts, simplifying the components, and normalizing their coefficients.
    %
    % Inputs:
    %   SymVarList - A list of symbolic variables (possibly complex).
    %
    % Outputs:
    %   SymVarList - The list of normalized symbolic variables.

    for i = 1:numel(SymVarList)
        % Extract real and imaginary parts of the symbolic variable
        A_SymVar = real(SymVarList(i));
        B_SymVar = imag(SymVarList(i));
        
        % Define a temporary symbolic variable for simplification
        syms zzzz real;
        
        % Get the children of the real and imaginary parts
        A_SymVar_children = children(A_SymVar + zzzz);
        B_SymVar_children = children(B_SymVar + zzzz);
        
        % Initialize final variables for real and imaginary parts
        A_SymVar_final = sym(0); 
        B_SymVar_final = sym(0);
        
        % Normalize the real part by processing its children
        for j = 1:length(A_SymVar_children) - 1
            A_SymVar_final = A_SymVar_final + iNormalizeSymVar(A_SymVar_children{j});
        end
        
        % Normalize the imaginary part by processing its children
        for j = 1:length(B_SymVar_children) - 1
            B_SymVar_final = B_SymVar_final + iNormalizeSymVar(B_SymVar_children{j});
        end
        
        % Recombine the normalized real and imaginary parts
        SymVarList(i) = A_SymVar_final + B_SymVar_final * 1i;
    end
end

function SymVar = iNormalizeSymVar(SymVar, Accuracy)
    % iNormalizeSymVar normalizes a symbolic variable by simplifying its
    % coefficients and matching them to known special values (such as integer multiples 
    % of square roots), returning the normalized expression.
    %
    % Inputs:
    %   SymVar   - The symbolic variable to normalize.
    %   Accuracy - A threshold to approximate coefficients (default is 1e-8).
    %
    % Outputs:
    %   SymVar   - The normalized symbolic variable.
    
    if nargin < 2
        Accuracy = 1e-8;  % Default accuracy if not provided
    end

    % Simplify the symbolic variable and check if it's zero
    if isequal(simplify(SymVar), sym(0))
        SymVar = sym(0);
        return;
    else
        % Extract the coefficients and base of the symbolic expression
        [Coeffs, SvmVarBase] = coeffs(SymVar);
        
        try
            % Attempt to convert coefficients to numeric values for comparison
            CoeffsNum = double(Coeffs);
        catch
            return;  % Return if conversion fails
        end
        
        % Known special bases and predefined fractions for coefficient comparison
        SpecialBase = [1, sqrt(2), sqrt(3), sqrt(5)];
        SpecialPre = [1, 2, 3, 4, 5, 6, 7, 8, 9, 1/2, 1/3, 1/4, 1/5, 2/3, 4/3, 8/3];
        
        % Compare the coefficients with special values
        for iSpecialBase = SpecialBase
            for jSpecialPre = SpecialPre
                CompareNum = jSpecialPre * iSpecialBase;
                
                % If coefficients are close to a special value, normalize accordingly
                if abs(CoeffsNum - CompareNum) < Accuracy
                    SymVar = sym(CompareNum) * SvmVarBase;
                    return;
                end
                if abs(CoeffsNum + CompareNum) < Accuracy
                    SymVar = -sym(CompareNum) * SvmVarBase;
                    return;
                end
            end
        end
    end
end
