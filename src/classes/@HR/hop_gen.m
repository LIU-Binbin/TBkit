function hop = hop_gen(hop_pre,nn_level)
% HOP_GEN Generate symbolic hopping term with level index
%
%   hop = HOP_GEN(hop_pre,nn_level) creates a symbolic hopping term
%   by appending the neighbor level to the base term.
%
%   INPUT ARGUMENTS:
%       hop_pre - Base hopping term (symbolic)
%       nn_level - Neighbor level index
%
%   OUTPUT ARGUMENTS:
%       hop - Symbolic hopping term with level suffix
%
%   SEE ALSO:
%       sym
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
tempstr = string(simplify(hop_pre));
tempstr = strcat(tempstr,"_",string(nn_level));
hop = sym(tempstr,'real');
end
