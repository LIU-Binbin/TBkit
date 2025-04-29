function C = minus(A,B)
%MINUS Subtraction operation for HK objects
%   C = MINUS(A, B) performs subtraction between HK objects or between
%   HK object and Term objects.
%
%   Inputs:
%       A - First operand (HK object or Term)
%       B - Second operand (HK object or Term)
%
%   Output:
%       C - Result of subtraction operation
%
%   See also HK, PLUS
if isa(A,'HK') && isa(B,'HK')
    C = A;
elseif isa(A,'HK') && ~isa(B,'HK')
    if isa(B,'Term')
        for i = 1:length(B)
            C = A.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
        end
    else
        C = A;
    end
elseif ~isa(A,'HK') && isa(B,'HK')
    if isa(A,'Term')
        if length(B) == 2
            for i = 1:length(A)
                C = B.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
            end
        else
            C = B;
        end
    else
        C = B;
    end
else
    C = 0;
end
end
