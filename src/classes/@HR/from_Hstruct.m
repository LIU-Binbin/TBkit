function H_hr = from_Hstruct(Hstruct)
% FROM_HSTRUCT Create HR object from structure array
%
%   H_hr = FROM_HSTRUCT(Hstruct) converts a structure array containing
%   Hamiltonian information into an HR object.
%
%   INPUT ARGUMENTS:
%       Hstruct - Structure array with fields:
%           vector: R-vectors (3-element vectors)
%           Hnum: Numerical Hamiltonian matrices
%           Hcoe: Symbolic Hamiltonian matrices
%           Degen: Degeneracy factors (optional)
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object constructed from input data
%
%   NOTES:
%       - Automatically determines WAN_NUM from first element
%       - Handles both numerical and symbolic Hamiltonian components
%       - Applies degeneracy factors if present
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

[NRPTS,~]=size(Hstruct);
try
WAN_NUM = length(Hstruct(1).Hnum);
catch
WAN_NUM = length(Hstruct(1).Hcoe);
end
V = [Hstruct.vector];
vectorL =(reshape(V,3,length(V)/3)');
HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
HcoeL = sym(zeros(WAN_NUM,WAN_NUM,NRPTS));
for i = 1:NRPTS
try
HnumL(:,:,i) = Hstruct(i).Hnum/Hstruct(i).Degen;
catch
HnumL(:,:,i) = zeros(WAN_NUM)  ;
end
try
HcoeL(:,:,i) = Hstruct(i).Hcoe;
catch
end
end
H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
end
