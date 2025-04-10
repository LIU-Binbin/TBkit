function Delta_Oper = Delta_Oper(Oper_list)
% DELTA_OPER Constructs a delta operator function from operator strings.
%
% SYNTAX:
%   Delta_Oper = Delta_Oper(Oper_list)
%
% DESCRIPTION:
%   This function processes a list of operator strings (Oper_list) to construct
%   a delta operator. Each operator string is expected to be in the format 
%   'Prefix_Direction_Offset', where Direction is 'x', 'y', or 'z' and Offset is 
%   either 'N' (for no offset) or a numerical value (as a string). The function 
%   builds a combined symbolic condition of the form:
%
%       (i1 == j1 + offset_x) & (i2 == j2 + offset_y) & (i3 == j3 + offset_z)
%
%   and converts it to a MATLAB function handle using matlabFunction. The resulting 
%   function has the signature @(i1,i2,i3,j1,j2,j3) and returns the evaluated delta condition.
%
% INPUT:
%   Oper_list - An array of operator strings, e.g., ["Oper_x_N", "Oper_y_2", "Oper_z_N"].
%
% OUTPUT:
%   Delta_Oper - A function handle that accepts integer indices (i1, i2, i3, j1, j2, j3)
%                and computes the delta condition based on the specified offsets.
%
% EXAMPLE:
%   Oper_list = ["Oper_x_N", "Oper_y_2", "Oper_z_N"];
%   Delta_Oper = Delta_Oper(Oper_list);
%   result = Delta_Oper(1, 2, 3, 1, 2, 3);
%
% SEE ALSO:
%   matlabFunction, sym, strsplit
%
% Author: [Your Name]
% Date: [Today's Date]
    tmpeq_all = "1";
    syms i1 i2 i3 j1 j2 j3 integer;
    tmpeqstr_list = ["(i1 == j1)", "(i2 == j2)", "(i3 == j3)"];
    
    for i = 1:length(Oper_list)
        tmpeq = "";
        strtmp = strsplit(string(Oper_list(i)), {'_'}, 'CollapseDelimiters', true);
        switch strtmp{2}
            case 'x'
                tmpeq = tmpeq + 'i1 == (j1';
            case 'y'
                tmpeq = tmpeq + 'i2 == (j2';
            case 'z'
                tmpeq = tmpeq + 'i3 == (j3';
            otherwise
        end
        if strcmp(strtmp{3}, 'N')
            tmpeq = tmpeq + ")";
        else
            tmpeq = tmpeq + " + " + strtmp{3} + ")";
        end
        tmpeq = "(" + tmpeq + ")";
        switch strtmp{2}
            case 'x'
                tmpeqstr_list(1) = tmpeq;
            case 'y'
                tmpeqstr_list(2) = tmpeq;
            case 'z'
                tmpeqstr_list(3) = tmpeq;
            otherwise
        end
    end
    
    tmpeq_all = tmpeq_all + "&" + tmpeqstr_list(1) + '&' + tmpeqstr_list(2) + '&' + tmpeqstr_list(3);
    Delta_Oper = matlabFunction(str2sym(tmpeq_all), 'Vars', [i1 i2 i3 j1 j2 j3]);
end
