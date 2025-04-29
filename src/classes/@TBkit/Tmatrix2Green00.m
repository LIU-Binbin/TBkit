function [Green_00,Green_00s1,Green_00s2]= Tmatrix2Green00(Tmatrix,H00,H01,w,eta,n)
% TMATRIX2GREEN00 Convert T-matrix to Green's functions for a 1D system
%
% Syntax:
%   [Green_00,Green_00s1,Green_00s2] = Tmatrix2Green00(Tmatrix,H00,H01,w,eta)
%   [Green_00,Green_00s1,Green_00s2] = Tmatrix2Green00(Tmatrix,H00,H01,w,eta,n)
%
% Description:
%   Computes the full, left surface, and right surface Green's functions
%   from a given T-matrix for a 1D system.
%
% Input Arguments:
%   Tmatrix - Input T-matrix
%   H00 - On-site Hamiltonian matrix
%   H01 - Nearest neighbor hopping matrix
%   w - Energy value (complex)
%   eta - Imaginary part for regularization
%   n - Power to raise T-matrix (default: 2)
%
% Output Arguments:
%   Green_00 - Full Green's function
%   Green_00s1 - Right surface Green's function
%   Green_00s2 - Left surface Green's function
%
% Example:
%   T = eye(2); H00 = diag([1,1]); H01 = eye(2);
%   [G,G1,G2] = Tmatrix2Green00(T,H00,H01,0.5,0.01);
if nargin < 6
    n=2;
end
Dimi = length(H00);
wc = w + 1i*eta;
%disp(Tmatrix);
Tmatrix = Tmatrix^n ;
eyemat = eye(length(Tmatrix));
%             if abs(det(Tmatrix)) <1e-6
%                 Tmatrix = Tmatrix + (1e-6)*eye(Dimi*2)*1i;
%             end
if rcond(Tmatrix) < 1e-10
    %[A,U] = eig(eyemat,Tmatrix);
    Tmatrix = Tmatrix + (1e-6)*eye(Dimi*2);
    [A,U] = eig(Tmatrix);
else
    [A,U] = eig(Tmatrix);
end
U_abs =abs(U);
[Asort,Usort] =park.sorteig(U_abs,A);
Lambda= diag(Usort);
S = Asort(:,Lambda<1);
SS = Asort(:,Lambda>1);
S2 = S(1:Dimi,:);
S1 = S(Dimi+1:2*Dimi,:);
S3 = SS(Dimi+1:2*Dimi,:);
S4 = SS(1:Dimi,:);
phi1 = H01*S2/S1;
phi2 = H01'*S3/S4;
Green_00 = inv(wc*eye(Dimi)-H00-phi1-phi2);
Green_00s2 = inv(wc*eye(Dimi)-H00-phi1);
Green_00s1 = inv(wc*eye(Dimi)-H00-phi2);
%     Green00.Green_00 =Green_00;
%     Green00.Green_00s1 =Green_00s1;
%     Green00.Green_00s2 =Green_00s2;
end