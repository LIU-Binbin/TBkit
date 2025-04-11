function H_hr = ForceTolist(H_hr)
% FORCETOLIST Convert HR object to list format
%
%   H_hr = FORCETOLIST(H_hr) converts HR object to list storage format.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object in list format
%
%   NOTES:
%       - Converts through full matrix format if sparse
%       - Uses rewrite method for final conversion
%
%   SEE ALSO:
%       HR, full, rewrite
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

switch H_hr.type
case 'sparse'
H_hr =H_hr.full();
H_hr = H_hr.rewrite;
case 'mat'
H_hr = H_hr.rewrite;
case 'list'
otherwise
error('Not support yet.');
end
end
