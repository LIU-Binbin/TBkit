function H_hr = ForceToMat(H_hr)
% FORCETOMAT Convert HR object to matrix format
%
%   H_hr = FORCETOMAT(H_hr) converts HR object to dense matrix format.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object in matrix format
%
%   NOTES:
%       - Converts sparse through full method
%       - Uses rewind for list format conversion
%
%   SEE ALSO:
%       HR, full, rewind
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

switch H_hr.type
case 'sparse'
H_hr =H_hr.full();
case 'mat'
case 'list'
H_hr = H_hr.rewind;
otherwise
error('Not support yet.');
end
end
