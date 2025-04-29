function C = lt(A,B)
%LT Less than comparison for HK objects
%   C = LT(A, B) performs element-wise comparison between HK objects or
%   between HK object and other types.
%
%   Inputs:
%       A - First operand (HK object or other)
%       B - Second operand (HK object or other)
%
%   Output:
%       C - Result of comparison operation
%
%   Note:
%       Currently supports limited operations between HK and char/double types
%
%   See also HK, LE
if isa(A,'HK') && isa(B,'HK')
    H_hk1 = A;
    H_hk2 = B;
    if H_hk1.Bassis_num ~= H_hk2.Bassis_num
        error('Bassis_num different');
    end
    error('not support at present.')
elseif isa(A,'HK') && ~isa(B,'HK')
    switch class(B)
        case 'char'
            switch B(1)
                case {'P','p'}
                    C = A.input_orb_struct(B,'vasp');
                case {'w','W'}
                case {'k','K'}
                    C = A.kpathgen3D(B);
                otherwise
            end
        case 'double'
            switch size(B,1)
                case A.WAN_NUN
                    C = A.input_orb_init(B);
                otherwise
            end
        otherwise
    end
elseif ~isa(A,'HK') && ~isa(B,'HK')
    error('not support at present.');
end
end
