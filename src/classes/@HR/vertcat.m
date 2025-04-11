function H_hr = vertcat(A,B)
%VERTCAT Vertical concatenation of HR objects
%
%   H_HR = VERTCAT(A, B) Vertically concatenates two HR objects
%   (actually performs horizontal concatenation via HORZCAT)
%
%   Inputs:
%       A, B - HR objects to concatenate
%
%   Output:
%       H_hr - Combined HR object
%
%   Note:
%       Currently just calls HORZCAT - may need proper vertical implementation
%
%   See also HORZCAT, CAT
H_hr = horzcat(A,B);
end
