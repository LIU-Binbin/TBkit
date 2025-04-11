% INIT Initialize the Htrig object with provided symbolic hopping terms.
%
% SYNTAX:
%   H_htrig = init(H_htrig, HsymL_trig_in, options)
%
% DESCRIPTION:
%   This function initializes a Htrig object by setting up its symbolic hopping term 
%   list (HsymL_trig) and corresponding coefficient matrix (HcoeL) based on the given 
%   input. The behavior depends on the object's Type:
%
%     - For 'sincos' type:
%         If HsymL_trig_in is not provided, a default list of symbolic expressions is used
%         (including sin and cos functions of k_x, k_y, and k_z). The function then generates 
%         the HcoeL field as a complex matrix (with real and imaginary parts given by symbolic 
%         variables 'A' and 'B'). The upper triangular part is symmetrized.
%
%     - For 'exp' type:
%         The function selects nearest-neighbor data from nn_store based on the specified 
%         level_cut option. If options.onsite_mode is set to 1, onsite hopping terms are generated 
%         for each basis element. Then, for each selected neighbor, hopping terms are constructed 
%         using symbolic variables (via SymbolicVarible) and an exponential factor expanded 
%         in k. Finally, the object is forced to be Hermitian.
%
% INPUTS:
%   H_htrig          - A Htrig object.
%   HsymL_trig_in    - (Optional) A symbolic array defining the hopping terms (default: 
%                      []; for 'sincos' type, a default set is used).
%   options          - A structure with the following optional fields:
%                         level_cut   - A double specifying the cutoff level for neighbor selection (default: 1).
%                         per_dir     - A vector specifying discretization per direction (default: [1,1,1]).
%                         onsite_mode - A flag (0 or 1) indicating whether to generate onsite hopping terms (default: 0).
%
% OUTPUT:
%   H_htrig          - The updated Htrig object with initialized HsymL_trig and HcoeL fields.
%
% EXAMPLE:
%   % Initialize a Htrig object H with default sincos type hopping terms:
%   H = init(H, [], struct('level_cut', 1, 'per_dir', [1,1,1], 'onsite_mode', 0));
%
function H_htrig = init(H_htrig, HsymL_trig_in, options)
arguments
    H_htrig Htrig;
    HsymL_trig_in sym = sym([]);
    options.level_cut double = 1;
    options.per_dir double = [1,1,1];
    options.onsite_mode double = 0;
end
syms k_x k_y k_z real;
k = [k_x; k_y; k_z];

if nargin < 2
    if strcmp(H_htrig.Type, 'sincos')
        HsymL_trig_in = [sym(1), sin(k_x), cos(k_x), sin(k_y), cos(k_y), sin(k_z), cos(k_z)];
    elseif strcmp(H_htrig.Type, 'exp')
        % For 'exp' type a default should be defined here if needed.
    end
end

if strcmp(H_htrig.Type, 'sincos')
    H_htrig.HsymL_trig = HsymL_trig_in;
    [sizeX, sizeY] = size(H_htrig.HcoeL);
    sizeH = [sizeX, sizeY, H_htrig.Kinds];
    HcoeL_tmp = sym('A', sizeH, 'real') + 1i*sym('B', sizeH, 'real');
    for i = 1:H_htrig.Kinds
        HcoeL_tmp(:, :, i) = triu(HcoeL_tmp(:, :, i)) + triu(HcoeL_tmp(:, :, i), 1)';
    end
    H_htrig.HcoeL = HcoeL_tmp;
elseif strcmp(H_htrig.Type, 'exp')
    select_nn_store = H_htrig.nn_store(H_htrig.nn_store(:,10) <= options.level_cut, :);
    nselect_nn_store = size(select_nn_store, 1);
    if options.onsite_mode == 1
        for i = 1:H_htrig.Basis_num
            A = sym(['A', '_', num2str(i), '_', num2str(i), '_0_ubar'], 'real');
            B = sym(['B', '_', num2str(i), '_', num2str(i), '_0_ubar'], 'real');
            SymHopping = A + 1i*B;
            SymVar = sym(1);
            H_htrig = H_htrig.set_hop(SymHopping, SymVar, [i, i]);
        end
    end
    for i = 1:nselect_nn_store
        A = Htrig.SymbolicVarible("A", select_nn_store(i, [6,7,8]), select_nn_store(i, 1:2), select_nn_store(i, 10));
        B = Htrig.SymbolicVarible("B", select_nn_store(i, [6,7,8]), select_nn_store(i, 1:2), select_nn_store(i, 10));
        SymHopping = A + 1i*B;
        SymVar = expand(exp(1i*select_nn_store(i, [3,4,5]) * k));
        H_htrig = H_htrig.set_hop(SymHopping, SymVar, select_nn_store(i, 1:2));
    end
    H_htrig = H_htrig.hermitize();
else
    % Additional types can be handled here.
end
end
