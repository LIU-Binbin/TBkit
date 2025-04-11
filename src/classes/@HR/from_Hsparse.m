function H_hr = from_Hsparse(Hsparse)
% FROM_HSPARSE Create HR object from sparse Hamiltonian
%
%   H_hr = FROM_HSPARSE(Hsparse) converts sparse Hamiltonian structure
%   to standard HR object.
%
%   INPUT ARGUMENTS:
%       Hsparse - Structure with sparse Hamiltonian data:
%           vectorL: R-vectors
%           HnumL: Cell array of sparse matrices
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object with dense matrix storage
%
%   NOTES:
%       - Converts all sparse matrices to full storage
%       - Creates empty symbolic component
%
%   SEE ALSO:
%       HR, full
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

vectorL = Hsparse.vectorL;
[NRPTS,~] = size(vectorL);
WAN_NUM = length(Hsparse.HnumL{1});
HcoeL = [];
HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
for i = 1:NRPTS
HnumL(:,:,i) = full(Hsparse.HnumL{i});
end
H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
end
