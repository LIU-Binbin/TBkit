function C = gt(B,A)
% GT Performs type-dependent operations for Hamiltonian objects or input configurations.
%   C = GT(B,A) handles operations between Htrig objects and other inputs (characters/numeric arrays). 
%   The function supports limited configurations and will error for most input combinations.
%
% Inputs:
%   B: 
%       - Htrig object: To be compared with A (only Basis_num check implemented).
%       - String character: Invokes specialized methods (e.g., 'P', 'K').
%       - Numeric array: Used for initialization when compatible with Htrig's properties.
%   A: 
%       - Htrig object: Central object for operations.
%       - Other types: Only compatible when B is also non-Htrig (currently unsupported).
%
% Outputs:
%   C: Result depends on input types:
%       - Symbolic orbital structure (when B is 'P'/'p').
%       - K-path generation structure (when B is 'k'/'K').
%       - Initialization structure (when B is numeric array matching A.WAN_NUN dimensions).
%       - Errors for unsupported cases.
%
% Supported Operations:
%   1. Htrig vs Htrig:
%      - Checks Basis_num compatibility (exact match required).
%      - Currently returns an error even when Basis_num matches.
%   2. Htrig vs character:
%      - 'P'/'p': Returns symbolic orbital representation
%      - 'k'/'K': Generates k-point path structure
%   3. Htrig vs numeric array:
%      - Requires array size to match A.WAN_NUN property
%      - Returns initialization-related structure
%
% Error Handling:
%   - Errors when:
%     * Basis_num mismatch between Htrig objects
%     * Unsupported input types or size mismatches
%     * Current implementation limitations
%   - Explicit warnings about incomplete functionality
%
% Notes:
%   - Non-standard (B,A) input order due to specific use-case requirements
%   - This function is under active development and has partial implementation
%   - Symbolic conversions applied to returned properties for 'P/p' cases
%   - Requires valid Htrig object properties for all operations (e.g., WAN_NUN, HsymL)
%
% Example Usage:
%   % Get symbolic orbital structure
%   C = gt('P', HtrigObj);
%   
%   % Generate k-path structure
%   C = gt('K', HtrigObj);
%   
%   % Initialize with parameter array (must match WAN_NUN dimensions)
%   params = [1 2 3];  % Example parameter array
%   C = gt(params, HtrigObj);
%
% See also: Htrig.input_orb_struct, Htrig.kpathgen3D

if isa(A,'Htrig') && isa(B,'Htrig')
    H_htrig1 = A;
    H_htrig2 = B;
    if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
        error('Bassis_num different');
    end
    error('not support at present.')
elseif isa(A,'Htrig') && ~isa(B,'Htrig')
    switch class(B)
        case 'char'
            switch B(1)
                case {'P','p'}
                    C = A.input_orb_struct(B,'sym');
                    C.Rm = sym(C.Rm );
                    C.orbL = sym(C.orbL );
                case {'w','W'}
                    % Currently unimplemented section
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
elseif ~isa(A,'Htrig') && ~isa(B,'Htrig')
    error('not support at present.');
end
end