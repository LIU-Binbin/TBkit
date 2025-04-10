% PLUS Overloads the addition operator (+) for Htrig objects.
%
% SYNTAX:
%   C = plus(A, B)
%
% DESCRIPTION:
%   This function overloads the addition operator for Htrig objects, allowing for 
%   the addition of two Htrig objects or the addition of a Htrig object with other 
%   types (e.g., Term, Trig, or symbolic expressions). The behavior depends on the 
%   types of the operands:
%
%   1. When both A and B are Htrig objects:
%         - The function extracts the HcoeL field from B, converts it into a cell array 
%           of matrices (one cell per "kind"), sets up the corresponding variable cells 
%           (default symbolic ones) and k-cell from B's HsymL_trig, and then calls the 
%           setup method on A to incorporate the data from B.
%
%   2. When A is a Htrig object and B is not a Htrig:
%         - If B is of class Term or Trig:
%               * The function updates the Htrig object A by adding B to its Trig_list,
%                 and then sets up rough symbolic-polynomial/hopping information using 
%                 B's properties.
%         - If B is symbolic (sym):
%               * The function treats B as a matrix of size [basis_num, basis_num] and 
%                 for each nonzero element, creates a temporary matrix with a single one 
%                 at the corresponding position. The symbolic element is then used to update 
%                 A via setup_rough.
%
%   3. When A is not a Htrig object and B is a Htrig object:
%         - The roles are reversed, and similar operations are performed with A's content.
%
%   4. In all other cases, the function returns one of the operands.
%
% INPUTS:
%   A - A Htrig object or a compatible type (Term, Trig, or sym).
%   B - A Htrig object or a compatible type (Term, Trig, or sym).
%
% OUTPUT:
%   C - The result of the addition operation, which is a Htrig object incorporating the 
%       contributions from A and B.
%
% EXAMPLE:
%   % Adding two Htrig objects:
%   C = H1 + H2;
%
%   % Adding a Htrig object with a symbolic matrix:
%   C = H + sym(someMatrix);
%
function C = plus(A, B)
    if isa(A, 'Htrig') && isa(B, 'Htrig')
        C = A;
        % Convert B.HcoeL into a cell array (one cell per kind)
        mat_cell = mat2cell(reshape(B.HcoeL, B.Basis_num, B.Basis_num * B.Kinds).', repmat(B.Basis_num, [1, B.Kinds]));
        Var_cell = mat2cell(sym(ones(1, B.Kinds)), 1);
        k_cell = mat2cell(B.HsymL_trig, 1);
        C = C.setup(Var_cell, k_cell, mat_cell);
        
    elseif isa(A, 'Htrig') && ~isa(B, 'Htrig')
        if isa(B, 'Term')
            C = A;
        elseif isa(B, 'Trig')
            C = A;
            C.Trig_list = C.Trig_list + B;
            for i = 1:length(B)
                C = C.setup_rough(B(i).symbolic_polynomial, B(i).pauli_mat);
            end
        elseif isa(B, 'sym')
            basis_num = A.Basis_num;
            if basis_num ~= length(B) || size(B,2) ~= size(B,2)
                error('Basis wrong!');
            end
            for i = 1:basis_num
                for j = 1:basis_num
                    if B(i,j) ~= sym(0)
                        tempmat = zeros(basis_num);
                        tempmat(i,j) = 1;
                        SymPoly = B(i,j);
                        A = A.setup_rough(SymPoly, tempmat);
                    end
                end
            end
            C = A;
        else
            C = A;
        end
        
    elseif ~isa(A, 'Htrig') && isa(B, 'Htrig')
        if isa(A, 'Term')
            C = B;
        elseif isa(A, 'Trig')
            C = B;
            C.Trig_list = C.Trig_list + A;
            for i = 1:length(A)
                C = C.setup_rough(A(i).symbolic_polynomial, A(i).pauli_mat);
            end
        else
            C = B;
        end
        
    else
        C = 0;
    end
end

