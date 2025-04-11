function H_hr = ForceToType(H_hr,Type)
% FORCETOTYPE Convert HR object to specified storage format
%
%   H_hr = FORCETOTYPE(H_hr,Type) converts HR object to requested format.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%       Type - Target type ('sparse','mat','list')
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object in requested format
%
%   NOTES:
%       - Wrapper for specific ForceTo* methods
%       - Provides unified interface for format conversion
%
%   SEE ALSO:
%       HR, ForceTosparse, ForceToMat, ForceTolist
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

switch Type
case 'sparse'
H_hr = H_hr.ForceTosparse();
case 'mat'
H_hr = H_hr.ForceToMat();
case 'list'
H_hr = H_hr.ForceTolist();
otherwise
error('Not support yet.');
end
end
