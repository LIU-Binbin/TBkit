function Gknew = CartisianMat(Gk,dir_seq,kstart)
%CARTISIANMAT Transform reciprocal lattice vectors to Cartesian basis
%
%   Syntax:
%       Gknew = CartisianMat(Gk,dir_seq,kstart)
%
%   Description:
%       Orthogonalizes reciprocal lattice vectors using Gram-Schmidt
%       process with specified direction priority.
%
%   Inputs:
%       Gk       - Input reciprocal lattice vectors [3×3]
%       dir_seq  - Direction priority sequence (default=[1,2,3])
%       kstart   - Starting direction ('kcar','k_x','k_y', or 'k_z')
%
%   Output:
%       Gknew    - Orthogonalized reciprocal lattice vectors
if nargin < 2
    dir_seq = [1,2,3];
end
if nargin < 3
    kstart = 'k_x';
end
switch kstart
    case 'kcar'
        Gknew = diag(diag(Gk));
    case 'k_x'
        k_x = Gk(dir_seq(1),:);
        Pk_x = TBkit.Pvector(k_x);
        if  dir_seq(2) == 0
            k_y = Gk(2,:);
        else
            k_y = Gk(dir_seq(2),:) - Gk(dir_seq(2),:)*Pk_x;
        end
        Pk_xy = TBkit.Pplane([k_x;k_y]);
        k_z = Gk(dir_seq(3),:) - Gk(dir_seq(3),:)*Pk_xy;
        Gknew =  [k_x;k_y;k_z];
    case 'k_y'
    case 'k_z'
        k_z = Gk(dir_seq(3),:);
        Pk_z = TBkit.Pvector(k_z);
        k_x = Gk(dir_seq(1),:) - Gk(dir_seq(1),:)*Pk_z;
        Pk_zx = TBkit.Pplane([k_z;k_x]);
        k_y = Gk(dir_seq(2),:) - Gk(dir_seq(2),:)*Pk_zx;
        Gknew =  [k_x;k_y;k_z];
        Gknew = Gknew(dir_seq,:);
end
% debug
%             V_Gk=dot(Gk(1,:),cross(Gk(2,:),Gk(3,:)));
%             V_Gknew=dot(Gknew(1,:),cross(Gknew(2,:),Gknew(3,:)));
%             disp([V_Gk,V_Gknew]);
end