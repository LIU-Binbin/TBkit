% LE Overloaded "less than or equal to" operator for Htrig objects.
%
% SYNTAX:
%   C = le(A, B)
%
% DESCRIPTION:
%   This function overloads the "less than or equal to" operator (<=) to handle
%   comparisons and operations involving Htrig objects. It supports the following cases:
%
%   1. Both A and B are Htrig objects:
%        - It first checks if their BASIS numbers (Bassis_num) are equal.
%        - If they differ, an error is thrown.
%        - If they are equal, the operation is not supported and an error is thrown.
%
%   2. A is a Htrig object and B is not:
%        - If B is of class char, the first character of B is inspected:
%              * If it is 'P' or 'p', the function calls input_orb_struct on A with B as parameter.
%              * If it is 'k' or 'K', the function calls kpathgen3D on A with B as argument.
%              * Cases for other characters (e.g., 'w', 'W') may be handled but are currently undefined.
%        - If B is of class double, the size of B in the first dimension is compared to A.WAN_NUN.
%              * If they match, the function calls input_orb_init on A with B.
%
%   3. Neither A nor B are Htrig objects:
%        - An error is thrown, as the operation is not supported.
%
% INPUTS:
%   A - A Htrig object or another type.
%   B - A Htrig object or another type (typically char or double when A is Htrig).
%
% OUTPUT:
%   C - The result of the operation; its meaning depends on the types of A and B.
%
% EXAMPLE:
%   % Example: Compare a Htrig object H with a character string:
%   C = H <= 'P_example';
%
function C = le(A, B)
    if isa(A, 'Htrig') && isa(B, 'Htrig')
        H_htrig1 = A;
        H_htrig2 = B;
        if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
            error('Bassis_num different');
        end
        error('not support at present.')
    elseif isa(A, 'Htrig') && ~isa(B, 'Htrig')
        switch class(B)
            case 'char'
                switch B(1)
                    case {'P', 'p'}
                        C = A.input_orb_struct(B, 'tbsym');
                    case {'w', 'W'}
                        % case for 'w' not implemented
                    case {'k', 'K'}
                        C = A.kpathgen3D(B);
                    otherwise
                        % Other cases can be added here
                end
            case 'double'
                switch size(B, 1)
                    case A.WAN_NUN
                        C = A.input_orb_init(B);
                    otherwise
                        % Other cases can be added here
                end
            otherwise
                % Other classes are not supported
        end
    elseif ~isa(A, 'Htrig') && ~isa(B, 'Htrig')
        error('not support at present.');
    end
end
