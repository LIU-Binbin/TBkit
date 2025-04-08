classdef Htrig < TBkit & matlab.mixin.CustomDisplay
properties
HcoeL;
HnumL;
HsymL_trig = sym(1);
end
properties
HsymL_coeL ;
HsymL_numL ;
Htrig_num   ;
Nslab = [0,0,0];
Htrig_latex ;
HsymL_trig_bk ;
seeds = [];
seedsvar = sym([]);
Sigmas        = []     ;
rm_list                ;
Trig_list =    Trig()  ;
Type {mustBeMember(Type,{'mat','list','sincos','exp','slab','sparse'})}= 'sincos';
Hmat_pre;
symvarL;
end
properties(Dependent = true)
Htrig_sym   ;
HsymL;
Kinds;
end
properties
num = false;
coe = true;
Duality_vector_dist;
end
methods (Access = protected)
 propgrp = getPropertyGroups(~)
end
methods
 H_htrig = Htrig(BASIS_NUM,Trig_list,options,propArgs)
end
methods
 H_htrig = setup_rough(H_htrig,symbolic_polynomial,pauli_mat,silence)
 H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence)
 H_htrig = set_hop(H_htrig,SymHopping,SymVar,Indij)
 [coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,symbolic_polynomial)
 H_htrig = find_HsymL_trig_bk(H_htrig,coeff_trig)
 Kind = k_symbol2Kind(H_htrig,k_symbol)
 H_htrig = reseq(H_htrig,wan_list,kinds_list)
 H_htrig = add_empty_one(H_htrig,vector)
end
methods
 H_htrig2 = rewrite(H_htrig,mode)
 H_htrig = diff(H_htrig,dir,options)
 H_htrig = simplify(H_htrig,Accuracy,options)
end
methods
 [HoutL] = HCAR_gen(H_htrig,klist,options)
 varargout = diff_klist(H_htrig,dir,klist,options)
end
methods
 H_htrig = init(H_htrig,HsymL_trig_in,options)
 H_htrig = applyOper(H_htrig,SymOper,options)
 H_htrig = dualize(H_htrig)
 H_htrig = hermitize(H_htrig)
 H_htrig_bk = subsOper(H_htrig,SymOper)
 H_htrig = applyU(H_htrig,U,conjugate ,antisymmetry )
 [H_htrig,H_htrig2] = applyR(H_htrig,R)
 [H_htrig,H_htrig2] = applyRt(H_htrig,Rt)
 [H_htrig,Smat] = Smatgen(H_htrig,R,Accuracy)
 Factorlist_parity = factorlist_parity(H_htrig)
 H_htrig = nn(H_htrig,search_range,Accuracy,Rlength_cut)
end
methods
 symvarL = get.symvarL(H_htrig)
 Htrig_sym = get.Htrig_sym(H_htrig)
 Htrig_latex = get.Htrig_latex(H_htrig)
 HsymL = get.HsymL(H_htrig)
 Kinds = get.Kinds(H_htrig)
end
methods
 C = plus(A,B)
 H_htrig = uminus(H_htrig)
 C = minus(A,B)
 C = mtimes(A,B)
 C = mtimes_inv(B,A)
 C = gt(B,A)
 C = lt(A,B)
 C = le(A,B)
 C = horzcat(A,B)
 H_htrig = ctranspose(H_htrig)
 H_htrig = conj(H_htrig)
 Htrig_sym = sym(H_htrig,options)
 Htrig_latex = latex(H_htrig)
end
methods
 H_htrig = translate(H_htrig,U)
 [H_htrig,EQL] = subsVar(H_htrig,varargin)
 H_htrig2 = subs(H_htrig,varargin)
 H_htrig = rotation(H_htrig,Rotation)
 H_htrig2 = discretize(H_htrig,Nslab,options)
 H_hr = Htrig2HR(H_htrig,options)
 H_hk = Htrig2HK(H_htrig,kpoints_f,options)
end
methods
 varargout = EIGENCAR_gen(H_htrig,options)
 [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_slab(H_htrig,options)
 [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_wire(H_htrig,options)
 H_htrig = Subsall(H_htrig,mode,AddtionVar)
 mat = HsymL_trig2mat(H_htrig,HsymL_trig)
 [num_label,coe_label,H_htrig] = NumOrCoe(H_htrig)
end
methods (Static)
 H_htrig = HR2Htrig(H_hr)
 Delta_Oper = Delta_Oper(Oper_list)
 sym_coe_list = coeff_extract(sym_term)
end
end
