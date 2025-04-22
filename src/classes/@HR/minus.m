function H_hr = minus(A,B)
%MINUS Subtraction operation for HR objects
%   Performs subtraction between HR objects or between HR and matrices.
%
%   Syntax:
%       H_hr = minus(A,B)
%
%   Inputs:
%       A - First operand (HR object or matrix)
%       B - Second operand (HR object or matrix)
%
%   Outputs:
%       H_hr - Resulting HR object from subtraction
%
%   Example:
%       H_diff = minus(H1, H2); % Subtract two HR objects
%       H_diff = minus(H1, eye(2)); % Subtract identity matrix from HR
H_hr = A + (-B);
end
