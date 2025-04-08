classdef HK < TBkit & matlab.mixin.CustomDisplay
    properties
        Degree = 1;
        Kinds;
        HcoeL;
        HnumL;
        HstrL;
        HsymL = sym([]);
        Hk_num   ;
    end
    properties
        HsymL_xyz = sym([]);
        nn_store_smart   ;
        nn_sparse_n      ;
        Atom_store_smart ;
        Rnn_map          ;
        Term_to_save     ;
        Trig_to_save     ;
        num = [];
        coe = [];
        Type = 'empty';
    end
    properties (Dependent=true,Transient = true)
        HsymL_k;
        Hk_sym ;
        Hk_latex ;
    end
    methods (Access = protected)
        propgrp = getPropertyGroups(~)
    end
    methods
        H_hk = HK(BASIS_NUM,Degree,Term_list,propArgs)

        
        H_hk = setup_single(H_hk,symbolic_polynomial,i,j,silence)
        H_hk = setup_rough(H_hk,symbolic_polynomial,pauli_mat,silence)
        H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence)
    end
    methods
        [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
        [H_sym_pauli,H_latex_pauli] = pauliDecomposition(H_hk)
    end
    methods
        HsymL_k = get.HsymL_k(H_hk)
        Hk_sym = get.Hk_sym(H_hk)
        Hk_latex = get.Hk_latex(H_hk)
        Type = get.Type(H_hk)
        [num_label,coe_label,H_hk] = NumOrCoe(H_hk)
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
    methods (Static,Hidden,Access= protected)
        str_out = quadk2K(str_in)
        [vector_list,Coeffs_list] = HstrL_classify(H_strL_tmp,R_struct,mode)
        [vector_list,Coeffs_list] = Hstr_mapping(H_str,R_struct,mode)
        [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
    end
end
