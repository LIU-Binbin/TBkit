function H_hr = sum(H_hr_list)
% SUM Sum multiple HR objects into one
%
%   H_HR = SUM(H_HR_LIST) combines multiple HR objects by summation
%
%   Input:
%       H_hr_list - Array of HR objects to sum
%   Output:
%       H_hr - Combined HR object
%
%   Notes:
%       - Uses overloaded plus operator
%       - Maintains all hopping terms
%       - Preserves storage format of first input
H_hr = H_hr_list(1);
for i = 2:length(H_hr_list)
H_hr = H_hr + H_hr_list(i);
end
end
