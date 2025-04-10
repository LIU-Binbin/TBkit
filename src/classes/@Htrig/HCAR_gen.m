function [HoutL] = HCAR_gen(H_htrig,klist,options)
% HCAR_GEN Generates Hamiltonian matrices for specified k-points based on system type.
%
% [HOUTL] = HCAR_GEN(H_htrig, klist, options) constructs a sequence of Hamiltonian matrices
% across a set of k-points, handling different Hamiltonian configurations (list/sparse/mat)
% and applying symmetry options as specified. The output format depends on the H_htrig.Type.
%
% Inputs:
%   H_htrig: Hamiltonian configuration structure/object with fields:
%       - Type: String defining the Hamiltonian representation ('list','sparse','mat').
%       - Basis_num: Integer number of basis functions in the system.
%       - HsymL_numL: Numeric matrix containing symmetry-related terms.
%       - HnumL: Cell/Matrix of numeric Hamiltonian contributions (varies by Type).
%       - Sparse_vector: Indices for sparse matrix storage (required for 'list'/'sparse' types).
%       - CutList: Bounds for summation ranges in sparse terms (required for 'list').
%       - HsymL_trig: Symbolic terms for sparse contributions (for 'sparse' type).
%       - Htrig_sym: Symbolic Hamiltonian expression (for default case).
%       - WAN_NUN: Property for numeric array initialization checks (observed in previous functions).
%   klist: Matrix of k-points (rows = individual k-points in Cartesian coordinates).
%   options: (Optional) Configuration structure with fields:
%       - Hermi (logical): Enforce Hermitian symmetry on outputs (default: true).
%
% Outputs:
%   HoutL: Hamiltonian matrices container depending on Type:
%       - 'list'/'mat': 3D numeric array (Basis_num x Basis_num x kn).
%       - 'sparse': Cell array of sparse matrices (each element corresponds to a k-point).
%       - Default: 3D numeric array with Hermitian symmetry enforced if options.Hermi is true.
%
% Notes:
%   - Requires valid H_htrig configuration for each Type (e.g., Sparse_vector for 'list' case).
%   - The 'Hermi' flag applies matrix symmetrization: (H + H')/2 for Hermitian matrices.
%   - Code contains experimental/unfinished sections (e.g., 'w/W' case in previous 'gt' function).
%   - Depends on external functions: SliceGen, TBkit.matrixtimespage.
%   - Potential code refinements:
%       * Check variable name spelling (e.g., 'Bassis_num' in earlier functions).
%       * Fix inconsistent assignments in non-Hermi branches (currently HoutL = HoutL).
%
% Example Usage:
%   % Generate Hermitian matrices for a 'list' type Hamiltonian
%   Hout = HCAR_gen(H_obj, k_points, struct('Hermi', true));
%   
%   % Sparse matrix generation with default options
%   H_out_cell = HCAR_gen(H_sparseObj, k_paths, struct());
%   
%   % Non-Hermitian mode for 'mat' type
%   Hout = HCAR_gen(H_matObj, k_list, struct('Hermi', false));
%
% See also: SliceGen, TBkit.matrixtimespage

arguments
    H_htrig;          % Hamiltonian configuration structure/object
    klist;            % k-points path (matrix with rows = Cartesian coordinates)
    options.Hermi = true;  % Symmetrize output matrices (Hermitian enforcement)
end

kn = size(klist,1);
switch H_htrig.Type
case 'list'
    if isempty(H_htrig.Sparse_vector) && isempty(H_htrig.CutList)
        H_htrig = SliceGen(H_htrig);
    end
    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
    Hnum_list_k = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:3)*klist.');
    for ki = 1:kn
        Hout = zeros(H_htrig.Basis_num,H_htrig.Basis_num);
        Hnum_list_ktmp = Hnum_list_k(:,ki);
        for i=1:H_htrig.N_Sparse_vector
            Hout(H_htrig.Sparse_vector(i,1),H_htrig.Sparse_vector(i,2)) = sum(Hnum_list_ktmp(H_htrig.CutList(i,1):H_htrig.CutList(i,2)));
        end
        if options.Hermi
            HoutL(:,:,ki) = (Hout+Hout')/2;
        else
            HoutL(:,:,ki) = HoutL;  % Potential assignment error (should be Hout?)
        end
    end
case 'sparse'
    for ki=1:kn
        Htemp=sparse(H_htrig.Basis_num,H_htrig.Basis_num);
        % Warning: ki variable potentially overwritten in inner loop
        for ki=1:H_htrig.Kinds
            Htemp = Htemp + Hnum_list{ki}*double(H_htrig.HsymL_trig(ki));
        end
        HoutL{ki} = Htemp;
    end
case 'mat'
    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
    for ki =1:kn
        Factorlist = exp(1i*H_htrig.HsymL_numL*klist(ki,:).');
        Hout = sum(TBkit.matrixtimespage(Factorlist,H_htrig.HnumL),3);
        if options.Hermi
            HoutL(:,:,ki) = (Hout+Hout')/2;
        else
            HoutL(:,:,ki) = HoutL;  % Potential assignment error (should be Hout?)
        end
    end
otherwise
    H_htrig.Htrig_num = subs(H_htrig.Htrig_sym);
    H_htrig.Hfun = matlabFunction(H_htrig.Htrig_num,'Vars',[sym('k_x'),sym('k_y'),sym('k_z')]);
    HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
    for ki =1:kn
        Hout = H_htrig.Hfun(klist(ki,1),klist(ki,2),klist(ki,3));
        if options.Hermi
            HoutL(:,:,ki) = (Hout+Hout')/2;
        else
            HoutL(:,:,ki) = HoutL;  % Potential assignment error (should be Hout?)
        end
    end
end
end