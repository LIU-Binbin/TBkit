function Equation_list = isolateAll(Equation_list,Symvar_list)
%ISOLATEALL Isolate variables in symbolic equations
%   EqList = isolateAll(EqList, SymList) attempts to isolate variables in
%   symbolic equations.
%
%   Inputs:
%       Equation_list - Array of symbolic equations
%       Symvar_list  - Target variables (default: variables in equations)
%
%   Output:
%       Equation_list - Modified equations with isolated variables
%
%   Handles different input dimensions and automatically simplifies equations
if nargin < 2
    Symvar_list = Equation_list;
end
zeroeq = sym(0) ==sym(0);
if isequal(size(Equation_list),size(Symvar_list))
    for i = 1:numel(Equation_list)
        Equation_list_tmp  = simplify(Equation_list(i));
        symvartmp = symvar(simplify(Symvar_list(i)));
        if  ~isequal(zeroeq,Equation_list_tmp) && ~isempty(symvartmp)
            try
                Equation_list(i) = isolate(Equation_list_tmp,symvartmp(1));
            catch
                fprintf('cant find the solution of Equation_list(%d):\n',i);
                disp(Equation_list(i));
            end
        end
    end
elseif isequal(size(Equation_list(:,:,1)),size(Symvar_list))
    %Nsym = numel(Symvar_list);
    SizeE = Equation_list ;
    for i = 1:numel(Equation_list)
        [row,col,~] = ind2sub(SizeE,i);
        symvartmp = symvar(Symvar_list(row,col));
        if Symvar_list(row,col)~=sym(0) && ~isequal(zeroeq,Equation_list(i))
            Equation_list(i) = isolate(Equation_list(i),symvartmp(1));
        end
    end
elseif isvector(Symvar_list)
    % not be implemented

else
    if Symvar_list~=sym(0)
        Symvar_list = symvar(Symvar_list);
        for i = 1:numel(Equation_list)
            if ~isequal(zeroeq,Equation_list(i))
                Equation_list(i) = isolate(Equation_list(i),Symvar_list);
            end
        end
    end

end
end

        