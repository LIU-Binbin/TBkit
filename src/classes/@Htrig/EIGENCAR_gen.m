% EIGENCAR_GEN Generate eigen energy and wavefunction carriers from a Htrig object.
%
% SYNTAX:
%   [EIGENCAR, WAVECAR] = EIGENCAR_GEN(H_htrig, options)
%
% DESCRIPTION:
%   This function computes the eigenvalues (EIGENCAR) and eigenfunctions (WAVECAR)
%   for a given Htrig object using different computational methods depending on the
%   type of the Htrig object (e.g., 'sparse', 'mat', 'list', 'sincos', etc.). 
%   It supports optional arguments to specify the Fermi energy, number of bands, 
%   k-point list, additional parameter substitution, display options, and progress 
%   printing.
%
%   For Htrig objects of type 'slab', the function prints a message and returns.
%
% INPUTS:
%   H_htrig         - A Htrig object.
%
%   options         - A structure with the following optional fields:
%       fermi       - Fermi energy (default: 0).
%       norb        - Number of enforced bands (default: -1, meaning use Basis_num).
%       klist       - k-point list (default: H_htrig.klist_cart).
%       para        - Parameter values for substitution (default: empty).
%       paraname    - Names of the parameters (no default).
%       show        - Logical flag to display the k-point path (default: false).
%       ax          - Handle for the axis for k-space plotting (default: []).
%       printmode   - Logical flag for enabling progress printing (default: true).
%
% OUTPUTS:
%   EIGENCAR        - Computed eigenenergies; for scalar parameter input, a matrix,
%                     or for multiple parameter sets, a cell array of matrices.
%   WAVECAR         - Computed eigenfunctions stored as wavefunction carriers.
%   (Optional) Additional output: The output from TBkit.klist_show if options.show is true.
%
% EXAMPLES:
%   % Example 1: Compute eigenenergies and eigenfunctions for a Htrig object H.
%   [EIGENCAR, WAVECAR] = EIGENCAR_GEN(H, struct('fermi', 0, 'norb', 20));
%
%   % Example 2: Compute and display eigenenergies along a specified k-path.
%   [EIGENCAR, WAVECAR, kplot] = EIGENCAR_GEN(H, struct('show', true, 'ax', gca));
%
% SEE ALSO:
%   TBkit.klist_show, eig, eigs, matlabFunction
%
function varargout = EIGENCAR_gen(H_htrig, options)
arguments
    H_htrig Htrig;
    options.fermi double = 0;
    options.norb double = -1;
    options.klist double = H_htrig.klist_cart;
    options.para  = [];
    options.paraname;
    options.show = false;
    options.ax = handle([]);
    options.printmode = true;
end

switch H_htrig.Type
    case 'slab'
        fprintf('use slab eigencar gen or ?\n');
        return;
    otherwise
end

fermi = options.fermi;
norb_enforce = options.norb;
if isempty(options.klist)
    H_htrig = H_htrig.kpathgen3D('KPOINTS');
    klist_cart_tmp = H_htrig.klist_cart;
else
    klist_cart_tmp = options.klist;
end

if options.show
    if isempty(options.ax)
        ax = TBkit.BZplot(H_htrig.Rm, 'color', 'r');
    else
        ax = options.ax;
    end
end

Hnum_list = H_htrig.HnumL;
if isempty(options.para)
    if H_htrig.Basis_num > 500
        print_mode = 1;
    else
        print_mode = 0;
    end
    [kn, ~] = size(klist_cart_tmp);
    if norb_enforce < 0
        NBANDS = H_htrig.Basis_num;
    elseif norb_enforce > 0
        NBANDS = norb_enforce;
    else
    end
    if strcmp(H_htrig.Type, 'list')
        H_htrig = H_htrig.SliceGen();
    end
    WAVECAR  = zeros(H_htrig.Basis_num, NBANDS, kn);
    EIGENCAR = zeros(NBANDS, kn);
    sizemesh = [H_htrig.Basis_num, H_htrig.Basis_num];
    if options.printmode
        pb = CmdLineProgressBar('BAND calculating ');
    end
    try
        HsymL_fun = matlabFunction(H_htrig.HsymL, 'Vars', H_htrig.VarsSeqLcart(1:H_htrig.Dim));
    catch
    end
    for ki = 1:kn
        switch H_htrig.Type
            case 'sparse'
                Htemp = sparse(H_htrig.Basis_num, H_htrig.Basis_num);
                for i = 1:H_htrig.Kinds
                    Htemp = Htemp + Hnum_list{i} * double(H_htrig.HsymL_trig(i));
                end
                Hout = Htemp;
            case {'mat'}
                Factorlist = exp(1i * H_htrig.HsymL_numL * klist_cart_tmp(ki, :).');
                Hout = sum(TBkit.matrixtimespage(Factorlist, H_htrig.HnumL), 3);
                Hout = (Hout + Hout') / 2;
            case 'list'
                Hout = zeros(sizemesh);
                Hnum_list_k = H_htrig.HnumL .* exp(1i * H_htrig.HsymL_numL(:, 1:H_htrig.Dim) * klist_cart_tmp(ki, :).');
                for i = 1:H_htrig.N_Sparse_vector
                    Hout(H_htrig.Sparse_vector(i, 1), H_htrig.Sparse_vector(i, 2)) = sum(Hnum_list_k(H_htrig.CutList(i, 1):H_htrig.CutList(i, 2)));
                end
                Hout = (Hout + Hout') / 2;
            case 'sincos'
                Input = num2cell(klist_cart_tmp(ki, :));
                kL = HsymL_fun(Input{:});
                Hout = H_htrig.HnumL;
                for i = 1:H_htrig.Kinds
                    Hout(:, :, i) = Hout(:, :, i) .* kL(:, i);
                end
                Hout = sum(Hout, 3);
                Hout = (Hout + Hout') / 2;
            otherwise
                Input = num3cell(klist_cart_tmp(ki, :));
                Hout = H_htrig.Hfun(Input{:});
                Hout = (Hout + Hout') / 2;
        end
        if norb_enforce < 0
            try
                [A, U] = eig(full(Hout));
            catch
                disp([ki, klist_cart_tmp(ki, :)]);
                disp(Hout);
                disp(H_htrig.Hfun);
                error('check this k point');
            end
        elseif norb_enforce > 0
            [A, U] = eigs(Hout, NBANDS, fermi);
            [A, U] = park.sorteig(U, A);
        else
        end
        EIGENCAR(:, ki) = diag(U);
        WAVECAR(:, :, ki) = A;
        if options.printmode
            pb.print(ki, kn, ' ...');
        end
    end
    if options.printmode
        pb.delete;
    end
else
    Npara = size(options.para, 1);
    paraN = size(options.para, 2);
    if H_htrig.Basis_num > 500
        print_mode = 1;
    else
        print_mode = 0;
    end
    [kn, ~] = size(klist_cart_tmp);
    if norb_enforce < 0
        NBANDS = H_htrig.Basis_num;
    elseif norb_enforce > 0
        NBANDS = norb_enforce;
    else
    end
    WAVECAR  = [];
    EIGENCAR{Npara} = zeros(NBANDS, kn);
    Htrig_sym = H_htrig.Hsym;
    for j = 1:Npara
        fprintf('**************************************************************************************\n');
        for i = 1:paraN
            fprintf('%s :', mat2str(string(sym(options.paraname(i)))));
            fprintf('%f\n', options.para(j, i));
        end
        EIGENCAR_tmp = zeros(NBANDS, kn);
        if strcmp(H_htrig.Type, 'sparse')
            for i = 1:H_htrig.Kinds
                Hnum_list{i} = subs(H_htrig.HnumL{i}, sym(options.paraname), options.para(j, :));
            end
        else
            H_fun_tmp = matlabFunction(subs(Htrig_sym, sym(options.paraname), options.para(j, :)), 'Vars', sym(["k_x", "k_y", "k_z"]));
        end
        for ki = 1:kn
            k_x = klist_cart_tmp(ki, 1);
            k_y = klist_cart_tmp(ki, 2);
            k_z = klist_cart_tmp(ki, 3);
            if strcmp(H_htrig.Type, 'sparse')
                Htemp = sparse(H_htrig.Basis_num, H_htrig.Basis_num);
                for i = 1:H_htrig.Kinds
                    Htemp = Htemp + Hnum_list{i} * double(H_htrig.HsymL_trig(i));
                end
                Hout = Htemp;
            else
                Hout = H_fun_tmp(k_x, k_y, k_z);
                Hout = (Hout + Hout') / 2;
            end
            if norb_enforce < 0
                try
                    [A, U] = eig(full(Hout));
                catch
                    disp([ki, k_x, k_y, k_z]);
                    disp(Hout);
                    disp(H_htrig.Hfun);
                    error('check this k point');
                end
            elseif norb_enforce > 0
                [A, U] = eigs(Hout, NBANDS, fermi);
                [A, U] = park.sorteig(U, A);
            else
            end
            EIGENCAR_tmp(:, ki) = diag(U);
            if print_mode == 1
                fprintf('%d th kpoints has been calculated in %d kpoints total\n', ki, kn);
            end
        end
        EIGENCAR{j} = EIGENCAR_tmp;
    end
end

varargout{1} = EIGENCAR;
varargout{2} = WAVECAR;
if options.show
    varargout{3} = TBkit.klist_show('klist', klist_cart_tmp, 'ax', ax);
end
end
