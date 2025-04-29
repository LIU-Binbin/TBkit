function A_1 = BerryConnection_definition(Eigenvector_sym,para1)
%BERRYCONNECTION_DEFINITION Calculate Berry connection for an eigenvector
%
%   Syntax:
%       A_1 = BerryConnection_definition(Eigenvector_sym,para1)
%
%   Description:
%       Computes Berry connection A = i<ψ|∂ψ/∂para> for a given
%       eigenvector and parameter.
%
%   Inputs:
%       Eigenvector_sym - Symbolic eigenvector
%       para1          - Parameter to differentiate with respect to
%
%   Output:
%       A_1 - Berry connection value
A_1 = 1i*Eigenvector_sym'*diff(Eigenvector_sym,para1);
end