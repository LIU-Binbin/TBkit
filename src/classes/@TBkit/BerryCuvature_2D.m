function [BC_2D,BCL] = BerryCuvature_2D(WAVECAR,sizemesh,options)
%BERRYCUVATURE_2D Calculate 2D Berry curvature from wavefunctions
%
%   Syntax:
%       [BC_2D,BCL] = BerryCuvature_2D(WAVECAR,sizemesh,options)
%
%   Description:
%       Computes Berry curvature on a 2D k-mesh using discrete wavefunction
%       overlaps. Implements the Fukui-Hatsugai-Suzuki method.
%
%   Inputs:
%       WAVECAR  - Wavefunction array (orbitals × bands × k-points)
%       sizemesh - Dimensions of k-mesh [Nk1,Nk2]
%       options  - Options structure:
%                  sum - Sum over bands flag (default=true)
%
%   Outputs:
%       BC_2D - Berry curvature array
%       BCL   - (Reserved for future use)
arguments
    WAVECAR
    sizemesh
    options.sum = true;
end
if options.sum
    Nband = 1;
else
    Nband = size(WAVECAR,2);
end
sizemesh_WAVECAR = sizemesh+1;% Contaion Edage % For intgrel method
kn = size(WAVECAR,3);
BC_2D = zeros([sizemesh,Nband]);
[i_list,j_list] = ind2sub(sizemesh_WAVECAR,1:kn);
% ind_VV_list = 1:kn;
% Vk
Vk_list = 1:kn;
seqL = i_list <=sizemesh(1) & j_list <=sizemesh(2);
Vk_list = Vk_list(seqL);
i_list = i_list(seqL);
j_list = j_list(seqL);
% Vk1
Vk1_i_list = i_list+1;%Vk1_i_list(Vk1_i_list>sizemesh(1)) = 1;
ind_Vk1_list = sub2ind(sizemesh_WAVECAR,Vk1_i_list,j_list);
% Vk2
Vk1_j_list = j_list+1;%Vk1_j_list(Vk1_j_list>sizemesh(2)) = 1;
ind_Vk2_list = sub2ind(sizemesh_WAVECAR,i_list,Vk1_j_list);
% Vk1k2
ind_Vk1k2_list = sub2ind(sizemesh_WAVECAR,Vk1_i_list,Vk1_j_list);
%
kn = length(Vk_list);
[iL,jL]= ind2sub(sizemesh,1:kn);
if options.sum
    for k = 1:kn
        VV = WAVECAR(:,:,Vk_list(k));
        Vk1 = WAVECAR(:,:,ind_Vk1_list(k));
        Vk2 = WAVECAR(:,:,ind_Vk2_list(k));
        Vk1k2 = WAVECAR(:,:,ind_Vk1k2_list(k));
        BC_2D(iL(k),jL(k),:) = TBkit.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
    end
else
    for i =1:Nband
        for k = 1:kn
            VV = WAVECAR(:,i,Vk_list(k));
            Vk1 = WAVECAR(:,i,ind_Vk1_list(k));
            Vk2 = WAVECAR(:,i,ind_Vk2_list(k));
            Vk1k2 = WAVECAR(:,i,ind_Vk1k2_list(k));
            BC_2D(iL(k),jL(k),i) = TBkit.BerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2);
        end
    end
end
if nargout == 2
    BCL = [];
end
end