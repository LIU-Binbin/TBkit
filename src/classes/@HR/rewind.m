function H_hr = rewind(H_hr)
% REWIND Convert HR object back to matrix storage format
%
%   H_HR = REWIND(H_HR) converts list-formatted HR object back to matrix format
%
%   Input:
%       H_hr - HR object in list format
%   Output:
%       H_hr - HR object in matrix format
%
%   Notes:
%       - Wrapper for rewrite() with rewind option
%       - Preserves all hopping information

H_hr = H_hr.rewrite('rewind',true);
end
