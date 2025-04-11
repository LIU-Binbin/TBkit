% LT Overloaded "less than" operator for Htrig objects.
%
% SYNTAX:
%   C = lt(A, B)
%
% DESCRIPTION:
%   This function overloads the "less than" operator (<) for operations involving 
%   Htrig objects. It handles several cases:
%
%   1. If both A and B are Htrig objects:
%          - The function checks that their BASIS numbers (Bassis_num) are equal.
%          - If they differ, an error is thrown.
%          - If they are equal, the operation is not yet supported and an error is raised.
%
%   2. If A is a Htrig object and B is not:
%          - When B is a char:
%                * If the first character is 'P' or 'p', the function calls input_orb_struct
%                  on A with the argument 'vasp'.
%                * If the first character is 'k' or 'K', the function calls kpathgen3D on A.
%                * Other cases (e.g., for 'w'/'W') may be added later.
%          - When B is a double:
%                * The function compares the number of rows of B with A.WAN_NUN. If they match,
%                  it calls input_orb_init on A with B.
%
%   3. If neither A nor B are Htrig objects:
%          - The operation is not supported and an error is thrown.
%
% INPUTS:
%   A - A Htrig object or another type.
%   B - A Htrig object or another type (typically char or double when A is Htrig).
%
% OUTPUT:
%   C - The result of the operation, whose meaning depends on the types of A and B.
%
% EXAMPLE:
%   % Compare a Htrig object H with a character string:
%   C = H < 'P_example';
%
function C = lt(A, B)
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
                    case {'P','p'}
                        C = A.input_orb_struct(B, 'vasp');
                    case {'w','W'}
                        % Operation for 'w' is not implemented.
                    case {'k','K'}
                        C = A.kpathgen3D(B);
                    otherwise
                        % Other cases may be added here.
                end
            case 'double'
                switch size(B,1)
                    case A.WAN_NUN
                        C = A.input_orb_init(B);
                    otherwise
                        % Not matching expected dimensions.
                end
            otherwise
                % Unsupported class for B.
        end
    elseif ~isa(A, 'Htrig') && ~isa(B, 'Htrig')
        error('not support at present.');
    end
end
