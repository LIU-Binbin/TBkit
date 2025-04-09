function C = lt(A,B)
%LT Overloaded less-than operator for HK objects
%
% Syntax:
%   C = A < B        % Operator form
%   C = lt(A,B)      % Functional form
%
% Inputs:
%   A - First HK object or comparison target
%   B - Second HK object or comparison specifier
%
% Output:
%   C - Result of comparison or operation result
%
% Description:
%   Handles special comparison cases for HK objects:
%   1. HK-HK comparison (currently unsupported)
%   2. HK-char comparison for path generation:
%      - 'P*'/'p*': Triggers input_orb_struct with 'vasp' mode
%      - 'K*'/'k*': Triggers kpathgen3D
%   3. HK-double comparison for orbital initialization
%
% Note:
%   Main difference from LE is default mode in input_orb_struct
%   Throws error for HK-HK comparisons
%
% Example:
%   Hk_new = Hk < 'P' % Equivalent to input_orb_struct('P','vasp')
%   Hk_path = Hk < 'k' % Generates k-path
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
