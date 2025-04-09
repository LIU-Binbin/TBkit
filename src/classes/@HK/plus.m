function C = plus(A, B)
%PLUS Overloaded addition for HK objects
%
% Syntax:
%   C = A + B       % Operator form
%   C = plus(A,B)   % Functional form
%
% Inputs:
%   A - HK object, Term, Trig, or symbolic matrix
%   B - HK object, Term, Trig, or symbolic matrix
%
% Output:
%   C - Combined Hamiltonian
%
% Description:
%   Handles multiple addition cases:
%   1. HK + HK: Coefficient-wise addition (requires matching basis/degree)
%   2. HK + Term: Adds terms using setup_rough
%   3. HK + Trig: Adds trigonometric terms
%   4. HK + sym: Adds symbolic matrix elements
%   5. Term + HK: Special handling for term addition
%
% Note:
%   Maintains Hamiltonian type (kp/tb) during addition
%   Enforces basis/dimension compatibility checks
%
% Example:
%   H_combined = Hk1 + Hk2; % Direct sum
%   H_with_term = Hk + Term(...); % Add new terms
if isa(A,'HK') && isa(B,'HK')
    H_hk1 = A;
    H_hk2 = B;
    if H_hk1.Basis_num ~= H_hk2.Basis_num
        error('basis num differ');
    end
    if H_hk1.Degree ~= H_hk2.Degree
        error('Degree differ');
    end
    C = H_hk1;
    C.HcoeL = C.HcoeL + H_hk2.HcoeL;
elseif isa(A,'HK') && ~isa(B,'HK')
    if isa(B,'Term')
        C = A;
        if contains(C.Type,'kp')
            C.Term_to_save = C.Term_to_save + B;
        else
            C.Term_to_save = B;
        end
        for i = 1:length(B)
            C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
        end
    elseif isa(B,'Trig')
        C = A;
        C.Trig_to_save = C.Trig_to_save+B.symbolic_polynomial*double(B.pauli_mat);
    elseif isa(B,'sym')
        basis_num = A.Basis_num;
        if basis_num~= length(B) || size(B,2)  ~=  size(B,2)
            error('Basis wrong!')
        end
        for i = 1:basis_num
            for j = 1: basis_num
                if B(i,j)~=sym(0)
                    tempmat = zeros( basis_num);
                    tempmat(i,j) =1 ;
                    A = A.setup_rough(B(i,j),tempmat);
                end
            end
        end
        C =A;
    else
        C = A;
    end
elseif ~isa(A,'HK') && isa(B,'HK')
    if isa(A,'Term')
        C = B;
        if contains(C.Type,'kp')
            C.Term_to_save = C.Term_to_save + A;
        else
            C.Term_to_save = A;
        end
        for i = 1:length(A)
            C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
        end
    elseif isa(A,'Trig')
        C = B;
        C.Trig_to_save = C.Trig_to_save + A.symbolic_polynomial*double(A.pauli_mat);
    else
        C = B;
    end
else
    C = 0;
end
end
