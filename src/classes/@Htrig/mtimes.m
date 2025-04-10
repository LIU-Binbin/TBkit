% MTIMES Overloads the multiplication operator (*) for Htrig objects.
%
% SYNTAX:
%   C = mtimes(A, B)
%
% DESCRIPTION:
%   This function defines the behavior of the multiplication operator (*) when one of 
%   the operands is a Htrig object. It supports scalar multiplication with either a 
%   double or symbolic (sym) scalar. The cases handled are:
%
%     1. Both A and B are Htrig objects:
%          - This operation is not implemented and an error is thrown.
%
%     2. A is a Htrig object and B is not:
%          - If B is a double, then for each "kind" in the Htrig object, HcoeL and HnumL 
%            are multiplied on the right by B.
%          - If B is symbolic, then HcoeL is multiplied on the right by B.
%
%     3. A is not a Htrig object and B is a Htrig object:
%          - If A is a double, then for each "kind", HcoeL and HnumL in B are multiplied 
%            on the left by A.
%          - If A is symbolic, then HcoeL is multiplied on the left by A.
%
%     4. In all other cases, the function returns 0.
%
% INPUTS:
%   A - Either a Htrig object or a scalar (double or sym).
%   B - Either a Htrig object or a scalar (double or sym).
%
% OUTPUT:
%   C - The product based on the type of A and B. For scalar multiplications, the internal 
%       fields (HcoeL and HnumL, if applicable) of the Htrig object are scaled accordingly.
%
% EXAMPLE:
%   % Multiply a Htrig object H by a double scalar:
%   C = H * 2;
%
%   % Multiply a symbolic scalar with a Htrig object H:
%   C = sym(3) * H;
%
% SEE ALSO:
%   mtimes, Htrig, sym, double
%
function C = mtimes(A, B)
    if isa(A, 'Htrig') && isa(B, 'Htrig')
        error('Multiplication of two Htrig objects is not implemented');
    elseif isa(A, 'Htrig') && ~isa(B, 'Htrig')
        C = A;
        if isa(B, 'double')
            for i = 1:C.Kinds
                C.HcoeL(:, :, i) = C.HcoeL(:, :, i) * B;
                C.HnumL(:, :, i) = C.HnumL(:, :, i) * B;
            end
        elseif isa(B, 'sym')
            for i = 1:C.Kinds
                C.HcoeL(:, :, i) = C.HcoeL(:, :, i) * B;
            end
        else
            % Do nothing or throw an error for unsupported classes.
        end
    elseif ~isa(A, 'Htrig') && isa(B, 'Htrig')
        C = B;
        if isa(A, 'double')
            for i = 1:C.Kinds
                C.HcoeL(:, :, i) = A * C.HcoeL(:, :, i);
                C.HnumL(:, :, i) = A * C.HnumL(:, :, i);
            end
        elseif isa(A, 'sym')
            for i = 1:C.Kinds
                C.HcoeL(:, :, i) = A * C.HcoeL(:, :, i);
            end
        else
            % Do nothing or throw an error for unsupported classes.
        end
    else
        C = 0;
    end
end
