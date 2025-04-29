function H_hk  = init(H_hk)
%INIT Initialize Hamiltonian coefficient matrix
%   H_hk = init(H_hk) creates a symbolic coefficient matrix for an HK object
%   with complex symmetric initialization.
%
%   Input/Output:
%       H_hk - HK object with initialized HcoeL property
%
%   The matrix is initialized as complex symmetric (A+iB where A is symmetric
%   and B is antisymmetric)
%
%   Example:
%       H = H.init(); % Initialize Hamiltonian coefficients
sizeH = size(H_hk.HcoeL);
HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
for i =1:H_hk.Kinds
HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
end
H_hk.HcoeL = HcoeL_tmp;
end
