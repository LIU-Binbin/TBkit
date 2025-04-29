function BC = BerryCuvature_Discrete_3D(VV,Vk1,Vk2,Vk3,Vk1k2,Vk2k3,Vk3k1)
%BERRYCUVATURE_DISCRETE_3D Calculate Berry curvature components in 3D
%
%   Syntax:
%       BC = BerryCuvature_Discrete_3D(VV,Vk1,Vk2,Vk3,Vk1k2,Vk2k3,Vk3k1)
%
%   Description:
%       Computes all three components of Berry curvature in 3D k-space
%       by calling BerryCuvature_Discrete_2D for each plane.
%
%   Inputs:
%       VV     - Wavefunction at base k-point
%       Vk1/2/3 - Wavefunctions at k-neighbors
%       Vk1k2 etc. - Wavefunctions at k-neighbor pairs
%
%   Output:
%       BC - Array of Berry curvature components [BC1,BC2,BC3]
BC(1) = TBkit.BerryCuvature_Discrete_2D(VV,Vk2,Vk3,Vk2k3);
BC(2) = TBkit.BerryCuvature_Discrete_2D(VV,Vk3,Vk1,Vk3k1);
BC(3) = TBkit.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
end