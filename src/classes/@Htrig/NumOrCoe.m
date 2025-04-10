% NUMORCOE Determine numerical and coefficient representation flags.
%
% SYNTAX:
%   [num_label, coe_label, H_htrig] = NumOrCoe(H_htrig)
%
% DESCRIPTION:
%   This function checks whether the Htrig object has its 'num' and 'coe'
%   flags defined. If either flag is empty, it calls the corresponding method 
%   in the TBkit superclass to determine these flags and updates the Htrig object.
%   Otherwise, it simply returns the existing flags.
%
% INPUT:
%   H_htrig  - A Htrig object.
%
% OUTPUT:
%   num_label  - A logical flag indicating whether the numerical representation is active.
%   coe_label  - A logical flag indicating whether the coefficient (symbolic) representation is active.
%   H_htrig    - The updated Htrig object with its 'num' and 'coe' fields set if they were empty.
%
% EXAMPLE:
%   [num_label, coe_label, H] = NumOrCoe(H);
%
function [num_label, coe_label, H_htrig] = NumOrCoe(H_htrig)
    if isempty(H_htrig.num) || isempty(H_htrig.coe)
        [num_label, coe_label] = NumOrCoe@TBkit(H_htrig);
        H_htrig.num = num_label;
        H_htrig.coe = coe_label;
    else
        num_label = H_htrig.num;
        coe_label = H_htrig.coe;
    end
end

