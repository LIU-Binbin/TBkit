% DIFF Compute the derivative of a Htrig object along a given direction.
%
%   H_htrig = DIFF(H_htrig, DIR, options) returns a Htrig object that is the 
%   derivative of the input Htrig object with respect to the variable associated 
%   with the index DIR. The method of differentiation depends on the Type of the 
%   Htrig object:
%
%     - 'sincos': Differentiates HsymL_trig symbolically, extracts and folds the 
%                 coefficient factors, and updates HnumL or HcoeL accordingly.
%     - 'exp':    Differentiates HsymL with respect to the k-space variable and 
%                 updates HnumL or HcoeL based on the ratio of the derivative to HsymL.
%     - 'mat':    Multiplies the appropriate column of HsymL_numL (or HsymL_coeL) by 1i 
%                 and applies it to HnumL (or HcoeL).
%     - 'list':   Performs element-wise 1i multiplication of HnumL or HcoeL with the 
%                 corresponding parts of HsymL_numL or HsymL_coeL.
%
%   Inputs:
%       H_htrig   - A Htrig object.
%       DIR       - (Optional) An integer specifying the direction for differentiation 
%                   (default: 1).
%       options   - (Optional) A structure with field:
%                      Accuracy - Tolerance for simplification (default: 1e-6).
%
%   Output:
%       H_htrig   - The updated Htrig object after differentiation.
%
%   Example:
%       % Differentiate a Htrig object along the first direction:
%       H = diff(H, 1, struct('Accuracy', 1e-6));
%
%   See also: simplify, matlab.sym, TBkit.matrixtimespage.
%
function H_htrig = diff(H_htrig,dir,options)
arguments
    H_htrig Htrig;
    dir = 1;
    options.Accuracy = 1e-6;
end

switch H_htrig.Type
    case 'sincos'
        VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
        K = VarUsing;
        HsymL_trig_tmp = diff(H_htrig.HsymL, K(dir));
        AdditionCoe = zeros(H_htrig.Kinds, 1, 'sym');
        for i = 1:H_htrig.Kinds
            [AdditionCoeTmp, HsymL_trigTmp] = coeffs(HsymL_trig_tmp(i));
            if isempty(AdditionCoeTmp) && isempty(HsymL_trigTmp)
                AdditionCoe(i) = 0;
                HsymL_trig_tmp(i) = 1;
            else
                AdditionCoe(i) = fold(@times, AdditionCoeTmp);
                HsymL_trig_tmp(i) = fold(@times, HsymL_trigTmp);
            end
        end
        if H_htrig.num
            H_htrig.HsymL_trig = HsymL_trig_tmp;
            H_htrig.HnumL = TBkit.matrixtimespage(double(AdditionCoe), H_htrig.HnumL);
        else
            H_htrig.HsymL_trig = HsymL_trig_tmp;
            H_htrig.HcoeL = TBkit.matrixtimespage(AdditionCoe, H_htrig.HcoeL);
        end
        
    case 'exp'
        syms k_x k_y k_z real;
        K = [k_x; k_y; k_z];
        AdditionCoe = diff(H_htrig.HsymL, K(dir)) ./ H_htrig.HsymL;
        if H_htrig.num
            H_htrig.HnumL = TBkit.matrixtimespage(double(AdditionCoe), H_htrig.HnumL);
        else
            H_htrig.HcoeL = TBkit.matrixtimespage(AdditionCoe, H_htrig.HcoeL);
        end
        
    case 'mat'
        if H_htrig.num
            H_htrig.HnumL = 1i * TBkit.matrixtimespage(H_htrig.HsymL_numL(:, dir), H_htrig.HnumL);
        else
            H_htrig.HcoeL = 1i * TBkit.matrixtimespage(H_htrig.HsymL_coeL(:, dir), H_htrig.HcoeL);
        end
        
    case 'list'
        if H_htrig.num
            H_htrig.HnumL = 1i * H_htrig.HnumL .* H_htrig.HsymL_numL(:, dir);
        else
            H_htrig.HcoeL = 1i * H_htrig.HcoeL .* H_htrig.HsymL_coeL(:, dir);
        end
end

H_htrig = simplify(H_htrig, options.Accuracy);
end
