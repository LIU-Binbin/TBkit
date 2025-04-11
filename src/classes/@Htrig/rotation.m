function H_htrig = rotation(H_htrig, Rotation)
%ROTATION Apply rotation transformation to the Htrig object.
%
%   H_htrig = ROTATION(H_htrig) rotates the symbolic momentum variables
%   of the Htrig object using the inverse of its reciprocal matrix Rm.
%
%   H_htrig = ROTATION(H_htrig, Rotation) applies a custom rotation matrix.
%   If Rotation is a string 'auto', the default rotation (inv(Rm)) is used.
%
%   Input:
%       H_htrig  - An object of class Htrig
%       Rotation - (optional) A symbolic rotation matrix (Dim x Dim) or the
%                  string 'auto' to trigger default behavior.
%
%   Output:
%       H_htrig  - Updated Htrig object with rotated symbolic variables.
%
%   Example:
%       H2 = rotation(H1);               % Use default inverse Rm rotation
%       H2 = rotation(H1, eye(3));       % Apply identity (no rotation)
%       H2 = rotation(H1, 'auto');       % Equivalent to default behavior
%
%   See also: subs

if nargin < 2
    Rotation = inv(sym(H_htrig.Rm));
end

if ischar(Rotation)
    if strcmp(Rotation,'auto')
        Rotation = inv(sym(H_htrig.Rm));
    end
end

VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
k = Rotation * VarUsing.';

for i = 1:H_htrig.Dim
    VarUsingCell{i} = VarUsing(i);
    kCell{i} = k(i);
end

H_htrig = H_htrig.subs(VarUsingCell, kCell);
end

