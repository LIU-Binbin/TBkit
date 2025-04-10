function vectorSeq = Getvector(H_hr,vector)
% GETVECTOR Find index of vector in HR object's vector list
%
%   vectorSeq = GETVECTOR(H_hr,vector) returns the index of the specified
%   vector in the HR object's vector list.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object containing vector list
%       vector - Vector to find (3-element array)
%
%   OUTPUT ARGUMENTS:
%       vectorSeq - Index of vector in H_hr.vectorL (0 if not found)
%
%   SEE ALSO:
%       HR, ismember
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
[~,vectorSeq]=ismember((vector),H_hr.vectorL,'rows');
end
