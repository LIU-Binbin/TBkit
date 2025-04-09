function H_hr = filter(H_hr,Accuracy)
% FILTER Remove small elements from Hamiltonian
%
%   H_hr = FILTER(H_hr,Accuracy) zeros out Hamiltonian elements smaller
%   than specified accuracy threshold.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format
%       Accuracy - Threshold for filtering (default: 1e-6)
%
%   OUTPUT ARGUMENTS:
%       H_hr - Filtered Hamiltonian
%
%   NOTES:
%       - Only works for 'list' and 'mat' Hamiltonian types
%       - Processes both real and imaginary parts separately
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
if nargin < 2
Accuracy = 1e-6;
end
switch H_hr.Type
case {'list','mat'}
HnumList = H_hr.HnumL;
HnumList_real = real(HnumList);
HnumList_imag = imag(HnumList);
HnumList_real(abs(HnumList_real) < Accuracy) = 0;
HnumList_imag(abs(HnumList_imag) < Accuracy) = 0;
H_hr.HnumL = HnumList_real + 1i*HnumList_imag;
otherwise
end
end
