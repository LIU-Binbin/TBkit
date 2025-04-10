% HERMITIZE Ensures that a Htrig object is Hermitian.
%
% SYNTAX:
%   H_htrig = hermitize(H_htrig)
%
% DESCRIPTION:
%   The hermitize function adjusts the internal matrix fields of a Htrig object so that
%   the object is Hermitian. For the coefficient representation (if coe_label is true),
%   it computes symbolic equations equating the real and imaginary parts of HcoeL with those
%   of its conjugate transpose, and substitutes these to enforce Hermiticity. For the numeric
%   representation (if num_label is true), it averages HnumL with its conjugate-transposed version.
%
% INPUT:
%   H_htrig - A Htrig object.
%
% OUTPUT:
%   H_htrig - The updated Htrig object with Hermitian matrix fields.
%
% EXAMPLE:
%   % Make the Htrig object H Hermitian:
%   H = hermitize(H);
%
function H_htrig = hermitize(H_htrig)
    [num_label, coe_label] = H_htrig.NumOrCoe();
    H_htrig_bk = H_htrig';
    if coe_label
        Equationlist_r = real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
        Equationlist_i = imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
        Equationlist_r = Htrig.isolateAll(Equationlist_r);
        Equationlist_i = Htrig.isolateAll(Equationlist_i);
        HcoeLtmp = subs(H_htrig.HcoeL, lhs(Equationlist_r), rhs(Equationlist_r));
        HcoeLtmp = subs(HcoeLtmp, lhs(Equationlist_i), rhs(Equationlist_i));
        H_htrig.HcoeL = HcoeLtmp;
    end
    if num_label
        H_htrig.HnumL = (H_htrig_bk.HnumL + H_htrig.HnumL) / 2;
    end
end
