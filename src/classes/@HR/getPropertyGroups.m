function propgrp = getPropertyGroups(~)
% GETPROPERTYGROUPS Define property groups for HR object display
%
%   propgrp = GETPROPERTYGROUPS(~) returns property groups for
%   customized display of HR objects.
%
%   OUTPUT ARGUMENTS:
%       propgrp - Property group object with standard properties
%
%   NOTES:
%       - Used by MATLAB's object display system
%       - Defines standard set of properties to show
%
%   SEE ALSO:
%       matlab.mixin.util.PropertyGroup
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

proplist = {'WAN_NUM','NRPTS','Type','HcoeL','HnumL','vectorL'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
