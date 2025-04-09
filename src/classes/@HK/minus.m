function C = minus(A,B)
%MINUS Overloaded subtraction for HK objects
%
% Syntax:
%   C = A - B       % Operator form
%   C = minus(A,B)  % Functional form
%
% Inputs:
%   A - HK object or Term array
%   B - HK object or Term array
%
% Output:
%   C - Result of operation
%
% Description:
%   Handles special cases of HK object subtraction:
%   1. HK - Term: Sets up terms using setup_rough
%   2. Term - HK: Special case for 2-element HK
%   3. Other combinations return unmodified input or 0
%
% Note:
%   Primary use is for Term object subtraction
%   Full HK-HK subtraction not currently implemented
%
% Example:
%   Hk_new = Hk - Term(...) % Adds negative terms
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
