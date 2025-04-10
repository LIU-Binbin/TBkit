function logicalal_num = ne(H_hr1,H_hr2)
% NE Overloaded not equal operator for HR objects
%
%   LOGICALAL_NUM = NE(H_HR1,H_HR2) compares two HR objects for inequality
%
%   Inputs:
%       H_hr1 - First HR object to compare
%       H_hr2 - Second HR object to compare
%   Output:
%       logicalal_num - Boolean result of comparison
%
%   Notes:
%       - Simply wraps the equality operator with logical negation
%       - Returns true if objects are not equal
logicalal_num = ~(H_hr1 == H_hr2);
end
