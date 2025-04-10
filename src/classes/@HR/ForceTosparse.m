function H_hr = ForceTosparse(H_hr)
% FORCETOSPARSE Convert HR object to sparse format
%
%   H_hr = FORCETOSPARSE(H_hr) converts HR object to sparse storage format.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object in sparse format
%
%   NOTES:
%       - Converts through matrix format if necessary
%       - Uses sparse method for final conversion
%
%   SEE ALSO:
%       HR, sparse
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
switch H_hr.type
case 'sparse'
case 'mat'
H_hr = H_hr.sparse;
case 'list'
H_hr = H_hr.rewind;
H_hr = H_hr.rewind;
otherwise
error('Not support yet.');
end
end
