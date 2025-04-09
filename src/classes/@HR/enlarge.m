function H_hr = enlarge(H_hr,dir,amp)
% ENLARGE Amplify Hamiltonian elements in specified direction
%
%   H_hr = ENLARGE(H_hr,dir,amp) multiplies Hamiltonian elements in the
%   specified direction by given amplitude factor.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format
%       dir - Direction ('x','y','z') or vector for selection
%       amp - Amplification factor
%
%   OUTPUT ARGUMENTS:
%       H_hr - Modified Hamiltonian with amplified elements
%
%   NOTES:
%       - When dir is char, selects elements based on orbital differences
%       - When dir is vector, selects elements directly
%       - Modifies the Hamiltonian in place
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

if isa(dir,'char')
H_hr = H_hr.rewrite;
vectorList = H_hr.vectorL;
orbList= H_hr.orbL;
orbListdiff = orbList(vectorList(:,H_hr.Dim+1),:) - orbList(vectorList(:,H_hr.Dim+2),:);
orbListdiff_r = orbListdiff*H_hr.Rm;
switch dir
case 'x'
selectL = logical(orbListdiff_r(:,1));
case 'y'
selectL = logical(orbListdiff_r(:,2));
case 'z'
selectL = logical(orbListdiff_r(:,3));
end
H_hr.HnumL(selectL)=H_hr.HnumL(selectL)*amp;
else
switch size(dir,2)
case 1
selectL = H_hr.vectorL(:,dir) ;
case 2
case 3
case 5
end
end
end
