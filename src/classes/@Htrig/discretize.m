% DISCRETIZE Convert a Htrig object to its slab representation.
%
% SYNTAX:
%   H_htrig2 = discretize(H_htrig, Nslab, options)
%
% DESCRIPTION:
%   The discretize method transforms a Htrig object into its slab representation. It 
%   optionally applies a rotation (specified by options.Rotation) and discretizes the 
%   system using the vector Nslab, which defines the number of layers (slab thickness) 
%   along each dimension. The function substitutes symbolic sine and cosine terms 
%   with corresponding sigma expressions, constructs the orbital lattice according 
%   to the discretization, and optionally removes orbitals via a user-provided removal 
%   function (options.rmfunc).
%
% INPUTS:
%   H_htrig   - A Htrig object.
%   Nslab     - A numeric vector defining the discretization of the slab (default: [0,10,0]).
%   options   - A structure with the following optional fields:
%                 rmfunc   - Function handle for orbital removal (default: @( ) 1).
%                 Rotation - Symbolic 3x3 rotation matrix (default: sym(eye(3))).
%
% OUTPUTS:
%   H_htrig2  - A new Htrig object in slab representation with updated properties,
%               including revised symbolic expressions (HsymL_trig), the orbital list (orbL),
%               and preliminary matrices (Hmat_pre).
%
% EXAMPLE:
%   % Convert a Htrig object H into a slab representation with a custom rotation and removal function:
%   opts.rmfunc = @myRemovalFunction;
%   opts.Rotation = sym([0 -1 0; 1 0 0; 0 0 1]);
%   H_slab = discretize(H, [0,12,0], opts);
%
% SEE ALSO:
%   rotation, subs, setup, fold, HsymL_trig2mat

function H_htrig2 = discretize(H_htrig, Nslab, options)
arguments
    H_htrig Htrig;
    Nslab double = [0,10,0];
    options.rmfunc function_handle = @()(1);
    options.Rotation  = sym(eye(3));
end

if strcmp(functions(options.rmfunc).function, '@()(1)')
    rm_mode = false;
else
    rm_mode = true;
end

if isequal(options.Rotation, sym(eye(3)))
    rotate_mode = false;
else
    rotate_mode = true;
end

if rotate_mode
    H_htrig = H_htrig.rotation(options.Rotation);
end

H_htrig2 = Htrig(H_htrig.Basis_num, 'Dim', H_htrig.Dim);
H_htrig2.seeds = ["Sigma_x", "Sigma_y", "Sigma_z", "Sigma_w"];
H_htrig2.Type = 'slab';
H_htrig2.Nslab = Nslab;
NSLAB = (H_htrig2.Nslab == 0) + H_htrig2.Nslab;
NS = fold(@times, NSLAB);
case_d = Nslab > 1;
seeds = ["x", "y", "z", "w"];
StrSinUsing = "sin(k_" + seeds + ")";
StrCosUsing = "cos(k_" + seeds + ")";
StrSigma_1__N  = "Sigma_" + seeds + "_1__N";
StrSigma_N__N  = "Sigma_" + seeds + "_N__N";
pat_d_pre   = "Sigma_" + seeds + "_";
SymSinUsing = str2sym(StrSinUsing);
SymCosUsing = str2sym(StrCosUsing);
SymSigma_1__N = str2sym(StrSigma_1__N);
SymSigma_N__N = str2sym(StrSigma_N__N);

for d = 1:numel(case_d)
    if case_d(d)
        H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig, SymSinUsing(d), 1i/2 * SymSigma_1__N(d));
        H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig, SymCosUsing(d), 1/2 * SymSigma_1__N(d));
        pat_d = pat_d_pre(d) + digitsPattern(1) + "__N";
        for i = 1:length(H_htrig.HsymL_trig)
            if ~contains(string(H_htrig.HsymL_trig(i)), pat_d)
                H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i) * SymSigma_N__N(d);
            end
        end
        H_htrig.Sigmas = [H_htrig.Sigmas; [SymSigma_N__N(d), SymSigma_1__N(d)]];
    end
end

H_htrig2.HsymL_trig_bk = [SymSigma_1__N SymSigma_N__N];
H_htrig2.HsymL_trig = sym([]);
count = 0;
for i = 1:H_htrig.Kinds
    [coeff_trig, symvar_list_trig, H_htrig2] = split_sym_eq(H_htrig2, H_htrig.HsymL_trig(i));
    for j = 1:numel(coeff_trig)
        count = count + 1;
        k_cell{count} = symvar_list_trig(j);
        mat_cell{count} = H_htrig.HcoeL(:, :, i);
        Var_cell{count} = coeff_trig(j);
    end
end
H_htrig2 = H_htrig2.setup(Var_cell, k_cell, mat_cell);
orb_tmp = zeros(NS * H_htrig.Basis_num, H_htrig.Dim);
NWAVE  = NS * H_htrig.Basis_num;
if isempty(H_htrig.orbL)
    H_htrig.orbL = zeros(H_htrig.Basis_num, H_htrig.Dim);
end
switch H_htrig.Dim
    case 1
        [i1L] = ind2sub(NSLAB, 1:NS);
        OrbAddL = i1L.' - 1;
    case 2
        [i1L, i2L] = ind2sub(NSLAB, 1:NS);
        OrbAddL = [i1L.' i2L.'] - 1;
    case 3
        [i1L, i2L, i3L] = ind2sub(NSLAB, 1:NS);
        OrbAddL = [i1L.' i2L.' i3L.'] - 1;
    case 4
        [i1L, i2L, i3L, i4L] = ind2sub(NSLAB, 1:NS);
        OrbAddL = [i1L.' i2L.' i3L.' i4L.'] - 1;
end
for i = 1:NS
    orb_tmp((i-1) * H_htrig.Basis_num + 1:i * H_htrig.Basis_num, :) = H_htrig.orbL + OrbAddL(i, :);
end
orb_tmp = orb_tmp ./ NSLAB;
SizeOrb_tmp = size(orb_tmp);
if rm_mode
    try
        for d = 1:SizeOrb_tmp(2)
            Input{d} = orb_tmp(:, d);
        end
        H_htrig2.rm_list = options.rmfunc(Input{:});
    catch
        H_htrig2.rm_list = false(1, NWAVE);
        for i = 1:NWAVE
            for d = 1:SizeOrb_tmp(2)
                Input{d} = orb_tmp(i, d);
            end
            H_htrig2.rm_list(i) = options.rmfunc(Input{:});
        end
    end
end
orb_tmp(H_htrig2.rm_list, :) = [];
H_htrig2.orbL = orb_tmp;
H_htrig2.Hmat_pre{numel(H_htrig2.HsymL_trig)} = sparse(NS, NS);
for i = 1:numel(H_htrig2.HsymL_trig)
    H_htrig2.Hmat_pre{i} = H_htrig2.HsymL_trig2mat(H_htrig2.HsymL_trig(i));
end
end
