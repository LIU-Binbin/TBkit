function A = GetProperty(H_hr,name)
% GETPROPERTY Get property value from HR object
%
%   A = GETPROPERTY(H_hr,name) returns the value of specified property.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to query
%       name - Property name to retrieve
%
%   OUTPUT ARGUMENTS:
%       A - Value of requested property
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

A = H_hr.(name);
end
