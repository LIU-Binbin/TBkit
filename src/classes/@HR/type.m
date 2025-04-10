function Type = type(H_hr)
%TYPE Get the storage type of HR object
%
%   TYPE = TYPE(H_HR) returns the storage format type of the HR object
%   ('mat', 'list', or 'sparse').
%
%   Input:
%       H_hr - HR object to query
%
%   Output:
%       Type - Storage type string
%
%   See also HR, CLASS
Type = H_hr.Type;
end
