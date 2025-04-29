function Bc = Berry_curvature_D2(Eigenvector_sym,para1,para2)
%BERRY_CURVATURE_D2 Calculate Berry curvature from eigenvector derivatives
%
%   Syntax:
%       Bc = Berry_curvature_D2(Eigenvector_sym,para1,para2)
%
%   Description:
%       Computes Berry curvature directly from eigenvector derivatives
%       using the formula: ∂A₂/∂para1 - ∂A₁/∂para2 where A is the
%       Berry connection.
%
%   Inputs:
%       Eigenvector_sym - Symbolic eigenvector
%       para1/2        - Parameters to differentiate with respect to
%
%   Output:
%       Bc - Berry curvature value
A_1 = Eigenvector_sym'*1i*diff(Eigenvector_sym,para1);
A_2= Eigenvector_sym'*1i*diff(Eigenvector_sym,para2);
Bc = simplify(diff(A_2,para1)-diff(A_1,para2));
end