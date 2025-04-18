classdef HK < TBkit & matlab.mixin.CustomDisplay
    %HK Class representing a Hamiltonian in k-space
    %
    %   This class handles creation and manipulation of k·p
    %   Hamiltonians in reciprocal space. It supports symbolic computation
    %   and various transformations.
    %
    % Properties:
    %   Degree - Order of k·p expansion
    %   Kinds - Number of unique terms in expansion
    %   HcoeL - Coefficient matrices for each term
    %   HnumL - Numeric coefficient matrices
    %   HstrL - String representations of terms
    %   HsymL - Symbolic representations of terms
    %   Hk_num - Numeric Hamiltonian
    %
    % Dependent Properties:
    %   HsymL_k - k-dependent symbolic terms
    %   Hk_sym - Symbolic Hamiltonian
    %   Hk_latex - LaTeX representation of Hamiltonian
    
    properties
        %DEGREE Order of k·p expansion
        %   Integer specifying maximum power of k in the Hamiltonian
        Degree = 1;
        
        %KINDS Number of unique terms in expansion
        %   Automatically determined from Degree
        Kinds;
        
        %HCOEL Coefficient matrices for each term
        %   3D array of symbolic matrices (basis×basis×Kinds)
        HcoeL;
        
        %HNUM Numeric coefficient matrices
        %   3D array of numeric matrices (basis×basis×Kinds)
        HnumL;
        
        %HSTRL String representations of terms
        %   Cell array of string representations
        HstrL;
        
        %HSYML Symbolic representations of terms
        %   Array of symbolic expressions for each term
        HsymL = sym([]);
        
        %HK_NUM Numeric Hamiltonian
        %   Numeric matrix representation of Hamiltonian
        Hk_num;
    end

    properties
        %HSYML_XYZ Symbolic terms in x,y,z coordinates
        %   Alternative symbolic representation
        HsymL_xyz = sym([]);
        
        %NN_STORE_SMART Storage for neighbor information
        nn_store_smart;
        
        %NN_SPARSE_N Sparse neighbor count
        nn_sparse_n;
        
        %ATOM_STORE_SMART Atomic information storage
        Atom_store_smart;
        
        %RNN_MAP Mapping of neighbor vectors
        Rnn_map;
        
        %TERM_TO_SAVE Saved symbolic terms
        Term_to_save;
        
        %TRIG_TO_SAVE Saved trigonometric terms
        Trig_to_save;
        
        %NUM Numerical parameters
        num = [];
        
        %COE Coefficient parameters
        coe = [];
        
        %TYPE Hamiltonian type identifier
        %   Can be 'empty', 'kp', 'tb', or 'kp&tb'
        Type = 'empty';
    end

    properties (Dependent=true,Transient = true)
        %HSYML_K k-dependent symbolic terms
        HsymL_k;
        
        %HK_SYM Symbolic Hamiltonian
        Hk_sym;
        
        %HK_LATEX LaTeX representation of Hamiltonian
        Hk_latex;
    end

    methods (Access = protected)
        %GETPROPERTYGROUPS Custom property display
        %   Overrides default property grouping for display
        propgrp = getPropertyGroups(~)
    end

    methods
        function H_hk = HK(BASIS_NUM,Degree,Term_list,propArgs)
            %HK Constructor for HK class
            %
            % Syntax:
            %   H_hk = HK()
            %   H_hk = HK(BASIS_NUM)
            %   H_hk = HK(BASIS_NUM,Degree)
            %   H_hk = HK(BASIS_NUM,Degree,Term_list)
            %   H_hk = HK(...,propArgs)
            %
            % Inputs:
            %   BASIS_NUM - Number of basis states (default=4)
            %   Degree - Order of k·p expansion (default=2)
            %   Term_list - Initial terms to add
            %   propArgs - Additional TBkit properties
            %
            % Output:
            %   H_hk - Constructed HK object
            arguments
                BASIS_NUM = 4;
                Degree = 2;
                Term_list = [];
                propArgs.?TBkit;
            end
            
            propArgsCell = namedargs2cell(propArgs);      
            H_hk = H_hk@TBkit(propArgsCell{:});
            switch nargin
                case 1
                    if isnumeric(BASIS_NUM)
                        H_hk.Degree = 1;
                        H_hk = H_hk.Degree2Kinds();
                        H_hk.Basis_num = BASIS_NUM;
                        H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Term_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Trig_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                    else
                        if isa(BASIS_NUM,'sym')
                            Hsym = BASIS_NUM;
                            if size(Hsym,1) ~= size(Hsym,2)
                                error('square sym mat required!')
                            end
                            try
                                H_hk.Degree = HK.checkDegree(Hsym);
                            catch
                                %H_hk.Degree = 1;
                            end
                            H_hk =H_hk.Degree2Kinds();
                            BASIS_NUM = length(Hsym);
                            H_hk.Basis_num = BASIS_NUM;
                            H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                            H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                            %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk.Term_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk.Trig_to_save =sym(zeros(BASIS_NUM,BASIS_NUM));
                            H_hk = H_hk + Hsym;
                        end
                    end
                case 2
                    if isnumeric(BASIS_NUM)
                        H_hk.Degree = Degree;
                        H_hk =H_hk.Degree2Kinds();
                        H_hk.Basis_num = BASIS_NUM;
                        H_hk.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        H_hk.HnumL = (zeros(BASIS_NUM,BASIS_NUM,H_hk.Kinds));
                        %H_hk.Hk_sym = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Term_to_save = sym(zeros(BASIS_NUM,BASIS_NUM));
                        H_hk.Trig_to_save = sym(zeros(BASIS_NUM,BASIS_NUM));
                    else

                    end
                case 3
                    if isnumeric(BASIS_NUM)
                        H_hk = HK(BASIS_NUM,Degree);
                        H_hk = H_hk + Term_list;
                    else

                    end

            end
        end

        
        H_hk = setup_single(H_hk,symbolic_polynomial,i,j,silence)
        H_hk = setup_rough(H_hk,symbolic_polynomial,pauli_mat,silence)
        H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence)
    end

    methods
        [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
        [H_sym_pauli,H_latex_pauli] = pauliDecomposition(H_hk)
    end

    methods
        function HsymL_k = get.HsymL_k(H_hk)
            %HSYML_K Get k-dependent symbolic terms
            %
            % Returns stored HsymL property (k-dependent terms)
            HsymL_k = H_hk.HsymL;
        end

        function Hk_sym = get.Hk_sym(H_hk)
            %HK_SYM Get symbolic Hamiltonian
            %
            % Constructs symbolic Hamiltonian from coefficient matrices and terms
            Hk_sym = sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
            for i =1:H_hk.Kinds
                Hk_sym = Hk_sym + H_hk.HcoeL(:,:,i)*H_hk.HsymL_k(i);
            end
            try
                Hk_sym = Hk_sym + H_hk.Trig_to_save; %temp
            catch
            end
        end

        function Hk_latex = get.Hk_latex(H_hk)
            %HK_LATEX Get LaTeX representation of Hamiltonian
            %
            % Converts symbolic Hamiltonian to LaTeX string
            Hk_latex = latex(H_hk.Hk_sym);
        end

        function Type = get.Type(H_hk)
            %TYPE Get Hamiltonian type identifier
            %
            % Determines if Hamiltonian is 'empty', 'kp', 'tb', or 'kp&tb'
            if isequal(H_hk.Trig_to_save,sym(zeros(size(H_hk.Trig_to_save))))
                TB = 0;
            else
                TB = 1;
            end
            if strcmp(class(H_hk.Term_to_save), class(sym(zeros(H_hk.Basis_num))))
                KP = 0 ;
            else
                KP =1;
            end
            if ~TB && ~KP
                Type = 'empty';
            elseif ~TB && KP
                Type = 'kp';
            elseif TB && ~KP
                Type = 'tb';
            elseif TB && KP
                Type = 'kp&tb';
            end
        end

        function [num_label,coe_label,H_hk] = NumOrCoe(H_hk)
            %NUMORCOE Get numerical and coefficient labels
            %
            % Returns stored numerical parameters and coefficients
            % If empty, inherits from vasplib superclass
            if isempty(H_hk.num) ||  isempty(H_hk.coe)
                [num_label,coe_label] = NumOrCoe@vasplib(H_hk);
                H_hk.num = num_label;
                H_hk.coe = coe_label;
            else
                num_label = H_hk.num;
                coe_label = H_hk.coe;
            end
        end
    end

    methods
        H_hk_out = PlaneWaveExpand(H_hk,N,ExpandDirection)
        H_hk = reseq(H_hk,basis_list)
    end

    methods
        C = plus(A,B)
        C = minus(A,B)
        C = mrdivide(A,B)
        C = mtimes(A,B)
        C = mtimes_inv(B,A)
        C = lt(A,B)
        C = le(A,B)
        C = horzcat(A,B)
        H_hk = conj(H_hk)
        H_hk = ctranspose(H_hk)
        H_hk_sym =sym(H_hk)
        [H_hk,Sublist,Unique_term] =  unique(H_hk,seed,checklist,options)
        H_hk = sum(H_hk_list)
        H_hk = simplify(H_hk,Accuracy)
    end

    methods
        H_hk = Degree2Kinds(H_hk)
    end

    methods
        H_hk  = init(H_hk)
        H_hk = applyU(H_hk,U,conjugate ,antisymmetry )
        matcell = matgen(H_hk,R,Accuracy)
        H_hk = applyR(H_hk,R)
        [H_hk_R,H_hk] = applyRU(H_hk,SymOper)
        [H_hk] = applyOper(H_hk,SymOper,options)
        H_hk = hermitize(H_hk)
        H_hk_bk = subsOper(H_hk,SymOper)
        J_mat = Jmat_gen(H_hk,QuantumL,options)
        H_hk = Subsall(H_hk,mode)
        H_hk  = subs(H_hk,varargin)
    end

    methods
        varargout= kp2TB(H_hk,kpoints_f,groups,options)
        varargout = EIGENCAR_gen(H_hk,options)
        Kind = k_symbol2Kind(H_hk,k_symbol)
    end

    methods(Static)
        Degree = checkDegree(Hsym,Dim)
        Orderlist  = orderlist(order)
        Factorlist_parity = factorlist_parity(Degree)
        sym_term = standardize_sym(sym_term)
        [sym_coe,sym_single,str_single] = coeff_extract(sym_term,Dim)
    end

    methods (Static,Hidden)
        str_out = quadk2K(str_in)
        [vector_list,Coeffs_list] = HstrL_classify(H_strL_tmp,R_struct,mode)
        [vector_list,Coeffs_list] = Hstr_mapping(H_str,R_struct,mode)
        [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
    end
end
