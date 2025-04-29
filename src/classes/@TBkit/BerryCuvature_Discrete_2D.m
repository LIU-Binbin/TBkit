function BC = BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2)
%BERRYCUVATURE_DISCRETE_2D Calculate Berry curvature for a single plaquette
%
%   Syntax:
%       BC = BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2)
%
%   Description:
%       Computes Berry curvature for a single k-space plaquette using
%       wavefunction overlaps. Based on the U(1) link variable approach.
%
%   Inputs:
%       VV    - Wavefunction at base k-point
%       Vk1   - Wavefunction at k1 neighbor
%       Vk2   - Wavefunction at k2 neighbor
%       Vk1k2 - Wavefunction at k1+k2 neighbor
%
%   Output:
%       BC - Berry curvature value (imaginary part of log of U(1) loops)

Uk1 = det(VV'*Vk1);Uk1 = Uk1/norm(Uk1);
Uk2 = det(VV'*Vk2);Uk2 = Uk2/norm(Uk2);
Uk1_k2 = det(Vk2'*Vk1k2);Uk1_k2 = Uk1_k2/norm(Uk1_k2);
Uk2_k1 = det(Vk1'*Vk1k2);Uk2_k1 = Uk2_k1/norm(Uk2_k1);
% ----- approach 2 -----
%             U = VV'*(Vk1*Vk1')*(Vk1k2*Vk1k2')*(Vk2*Vk2')*VV;
%             BC = angle(eig(U));
% berry curvature
BC = log(Uk1*Uk2_k1/(Uk1_k2*Uk2));
BC = imag(BC);
%BC = sum(BC);
end