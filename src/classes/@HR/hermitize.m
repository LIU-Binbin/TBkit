function H_hr = hermitize(H_hr)
% HERMITIZE Enforce Hermiticity on HR object
%
%   H_hr = HERMITIZE(H_hr) ensures the Hamiltonian is Hermitian by
%   averaging with its conjugate transpose.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to hermitize
%
%   OUTPUT ARGUMENTS:
%       H_hr - Hermitized HR object
%
%   NOTES:
%       - Handles both symbolic and numerical components
%       - Preserves vector hopping representation
%       - Uses symbolic solving for coefficient matching
%
%   SEE ALSO:
%       HR, isolateAll
%
H_hr_bk = H_hr';
if H_hr.vectorhopping
    H_hr = H_hr+H_hr_bk;
    H_hr = H_hr.simplify;
else
    if H_hr.coe
        Equationlist_r = real(H_hr.HcoeL - H_hr_bk.HcoeL) == 0;
        Equationlist_i = imag(H_hr.HcoeL - H_hr_bk.HcoeL) == 0;
        Equationlist_r = HR.isolateAll(Equationlist_r);
        Equationlist_i = HR.isolateAll(Equationlist_i);
        HcoeLtmp = subs(H_hr.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
        HcoeLtmp = subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));
        H_hr.HcoeL = HcoeLtmp;
    end
    if H_hr.num
        H_hr.HnumL = (H_hr_bk.HnumL + H_hr.HnumL )/2;
    end
end
end
