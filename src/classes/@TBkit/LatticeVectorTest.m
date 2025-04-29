function TrueOrFalse = LatticeVectorTest(Vector,Accuracy)
%LATTICEVECTORTEST Check if vector components are integers within tolerance
%   TF = LATTICEVECTORTEST(Vector, Accuracy) tests if all components of the
%   input vector are integers within specified accuracy.
%
%   Inputs:
%       Vector    - Input vector to test
%       Accuracy  - Tolerance for integer check (default: 1e-6)
%
%   Output:
%       TrueOrFalse - Logical true if all components are integers within tolerance
%
%   Example:
%       tf = LatticeVectorTest([1.000001, 2.0, -3.000003], 1e-5)
if nargin < 2
    Accuracy = 1e-6;
end
TrueOrFalse = true;
for i = 1:numel(Vector)
    if abs(rem(Vector(i),1)) > Accuracy
        TrueOrFalse = false;
        return;
    end
end
end