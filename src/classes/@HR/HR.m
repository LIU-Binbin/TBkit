classdef HR <TBkit & matlab.mixin.CustomDisplay
properties
vectorL ;
HnumL   ;
HcoeL   ;
end
properties
Duality_vector_dist;
end
properties
Type;
overlap logical= false;
num logical= false;
coe logical= true;
soc logical= false;
AvectorL;
BvectorL;
CvectorL;
vectorhopping = false;
end
properties (Transient,Hidden = true)
R_vector_dist;
end
properties(Dependent = true)
NRPTS    ;
WAN_NUM  ;
Line_000 ;
homecell ;
end
properties(Dependent = true,Hidden = true)
end
properties (Hidden = true)
nn_store_smart   ;
nn_sparse_n      ;
Atom_store_smart ;
Rnn_map          ;
ScoeL            ;
SnumL            ;
vectorL_overlap  ;
end
methods (Access = protected)
 propgrp = getPropertyGroups(~)
end
methods
 H_hr = HR(WAN_NUM,vectorL,options,propArgs)
end
methods (Static)
 H_hr = from_Hstruct(Hstruct)
 H_hr = from_Hsparse(Hsparse)
 H_hr = from_POSCAR_SE(POSCAR_file,options)
 H_hr = from_hdf5(filename)
 H_hr = from_wannier90(filename,Type,options)
end
methods(Static,Hidden,Access= protected)
 [dataArray,NRPT_list,NRPTS,NUM_WAN]=hrdat_read(filename)
end
methods
 H_hr = add_empty_one(H_hr,vector)
 H_hr = expand_empty_one(H_hr,orbOne,QuantumOne,elementOne)
 H_hr = set_hop(H_hr,amp,hi,hj,vector_list,mode)
 H_hr = set_hop_mat(H_hr,amp,vector,mode)
 H_hr = set_hop_single(H_hr,amp,hi,hj,vector,mode)
 H_hr = set_overlap(H_hr,amp,si,sj,vector_list,mode)
 H_hr = set_overlap_mat(H_hr,amp,vector,mode)
 H_hr = set_overlap_single(H_hr,amp,si,sj,vector,mode)
end
methods
 H_hr = set.NRPTS( H_hr, ~ )
 H_hr = set.WAN_NUM(H_hr,~)
end
methods
 NRPTS = get.NRPTS(H_hr)
 WAN_NUM = get.WAN_NUM(H_hr)
 Line_000 = get.Line_000(H_hr)
 homecell = get.homecell(H_hr)
 Type = type(H_hr)
end
methods
 vectorSeq = Getvector(H_hr,vector)
 H_hr = autohermi(H_hr,mode,options)
end
methods
 H_hr = plus(A,B)
 H_hr = minus(A,B)
 H_hr = uminus(H_hr)
 H_hr = times(A,B)
 H_hr = mrdivide(A,B)
 H_hr = premtimes(A,B)
 H_hr = mtimes(A,B)
 H_hr = power(H_hr,b)
 H_hr = mpower(A,B)
 varargout = eq(H_hr1,H_hr2)
 logicalal_num = ne(H_hr1,H_hr2)
 H_hr = ctranspose(H_hr)
 H_hr = transpose(H_hr)
 H_hr = horzcat(A,B)
 H_hr = vertcat(A,B)
 C = gt(B,A)
 C = lt(A,B)
 C = le(A,B)
 H_hr = kron(A,B)
 H_hr = conj(H_hr)
 H_hr = full(H_hr)
 H_hr = sparse(H_hr)
 [H_hr,EQL] = subs(H_hr,varargin)
 H_hr = simplify(H_hr,Accuracy)
 H_hr = filter(H_hr,Accuracy)
 Hsym = sym(H_hr,options)
 H_hr = sum(H_hr_list)
 [H_hr,Sublist,Unique_term] = unique(H_hr,seed,checklist,Accuracy)
end
methods
H_hr = Hnanowire_gen(H_hr,Nslab,np,vacuum_mode,options);
 H_hr = reseq(H_hr,wan_list,nrpt_list,nrpt_list_S)
 H_hr = cut_orb(H_hr,rm_list,options)
 H_hr = clean(H_hr,WANNUM)
 H_hr = project(H_hr,BASIS_MAT)
 H_hr = charalize(H_hr)
 H_hr = ForceToMat(H_hr)
 H_hr = ForceTosparse(H_hr)
 H_hr = ForceTolist(H_hr)
 H_hr = ForceToType(H_hr,Type)
 H_hr = OpenBoundary(H_hr,OBC_list)
end
methods
 H_hr = enlarge(H_hr,dir,amp)
 H_hr = add_soc(H_hr)
 H_hr = addorb(H_hr,orblist,options)
 H_hr = add_orb(H_hr,hop_struct,orbOne,QuantumOne,elementOne)
 H_hr = addsoc(H_hr,quantumL)
 H_hr = deltarule(H_hr,level_cut,mode,options)
 H_hr = alpharule(H_hr,level_cut,mode,options)
end
methods
 H_hr = cut_piece(H_hr,repeatnum,fin_dir,glue_edges,vacuum_mode)
 H_hr = supercell_hr(H_hr,Ns,options)
 H_hr = unfold_hr(H_hr,Ns,options)
 H_hr = translation(H_hr,translation_vector,options)
 [sc_orb,sc_vec,sc_elementL,sc_quantumL] = supercell_orb(H_hr,Ns,Accuracy)
 [pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,orb_id_L,pc_orb_id_L,pc_orb_selectL] = unfold_orb(H_hr,Ns,Accuracy,orb_id_L)
 H_hr = descritize(H_hr,Nslab,options)
 orbital_out = nanowire_orb(H_hr,fin_dir,vacuum_mode,options)
 H_hr = supercell(H_hr,Ns,filename,Rm,sites,Atom_name,Atom_num,findir)
 [H_hr_out,H_hr_pi_plus,H_hr_pi_minus] = realmap(H_hr)
end
methods (Static,Hidden,Access= protected)
 Hnum_list_wire_iz =  Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab,Type)
 Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
 cutlist= unique_label2cutlist(unique_label,NRPTS)
end
methods
 H_hr= hrz_gen(H_hr,kpoints_f,fin_dir,mode)
 [H00_H11_cell_list_1,H00_H11_cell_list_2] = H00_H11_cell_list_gen(H_hr,fin_dir,principle_layer)
 varargout = Green_prepare(H_hrz,principle_layer,fin_dir)
 [klist_cart,klist_frac,gap_list,fig] = findnodes2(H_hr,kzf_list,Noccupy,tolerance)
 H_hr = Subsall(H_hr,mode)
 H_hr = input_orb_struct(H_hr,filename,mode,options)
end
methods
 H_hr = GenfromOrth(H_hr,seed_r,seed_i,Accuracy,options)
 hrdat = Gen_hr(H_hr,filename,mode)
 [H_hr_bk] = POSCAR_gen(H_hr,filename,options)
 H_hr = tbbox_in_gen(H_hr,options)
 H_hr = pythtb_gen(H_hr,filename,options)
end
methods
 HcktObj = HR2Hckt(H_hr,options,options_homecell,options_para)
 [H_hr_forHckt,maxOnsite] = HRforHckt(H_hr,options)
 H_htrig = HR2Htrig(H_hr,options)
 H_hk = HR2HK(H_hr,kpoints_frac,options)
end
methods
varargout = EIGENCAR_gen(H_hr,options)
EIGENCARout = EIGENCAR_gen_sparse(H_hr,fermi,norb_enforce,klist_s_tmp)
end
methods
 H_hr = rewrite(H_hr,options)
 H_hr = rewind(H_hr)
 H_hr = init(H_hr,options)
 H_hr = applyOper(H_hr,SymOper,options)
 H_hr = symmetrize(H_hr,SymOper,options)
 H_hr = dualize(H_hr)
 [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf)
 [ml_cell,ij_list] = Tij2lm(H_hr,Rf)
 [H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper)
 H_hr = hermitize(H_hr)
 H_hr_bk = subsOper(H_hr,SymOper)
 [H_hr_R,H_hr] = applyRU(H_hr,SymOper )
 [H_hr_R,H_hr] = applyR(H_hr,R)
 [H_hr,Smat] = Smatgen(H_hr,R,Accuracy)
 Factorlist_parity = factorlist_parity(H_hr)
 H_hr = nn(H_hr,search_range,Accuracy,Rlength_cut,options)
end
methods
 [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_gen_3D(H_hr,kmesh,k3d,options)
 [EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = slab(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
 [EIGENCAR_slab,WEIGHTCAR_slab,klist1,klist2] = slab_fermiarc(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
 [EIGENCAR,orb_list,WEIGHTCAR,klist_l,kpoints_l,kpoints_name] = EIGENCAR_gen_wire(H_hr,Nslab,fermi,norb_enforce,KPOINTS_wire,vacuum_mode,np)
 [EIGENCAR,orb_list,WAVECAR] = EIGENCAR_gen_disk(H_hr,Nslab,fermi,norb_enforce,kpoints,vacuum_mode,np)
 [DOSCAR_l,DOSCAR_b,DOSCAR_r,w_list,klist_l,kpoints_l,kpoints_name] = surf(H_hr,w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode)
 [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2] = fermiarc(H_hr,w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode)
 [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2,E_list] = fermiarc3D(H_hr,w_range,fin_dir,kmesh,kfermiarc,options)
end
methods(Static)
 H_hr = from_POSCAR_Script(varargin)
end
methods
 A = list(H_hr,vectorL,options)
 H_hr = rmhop(H_hr,vectorL,options)
 [fig,ax] = show(H_hr,mode,options)
 Hout = printout(H_hr,print_list,mode)
 [VarInit,EQL2] =GetInit(H_hr,H_hr2,vectorL)
 H_atom_soc = H_atom_soc(H_hr)
 kloop = kloop_gen(H_hr,input_mat,mode)
 stucture
 [Rnn,nn_store_smart,Atom_store_smart,Rnn_map] = nn_information(H_hr,silence)
end
methods
 H_hr = H_TBSK_gen(H_hr,options)
 H_hr = H_TBSK_gen_sparse(H_hr,options)
end
methods
 H_hr_n = Connect2D(H_hr_n,H_hr_unitcell,opt)
end
methods(Static,Hidden,Access= protected)
 varargout = TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,options)
 varargout = TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,options)
 varargout = TBSK_Var_gen_single(L_1,L_2,m_1,m_2,nn_level,l,m,n,options)
 Coeff = TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,options)
end
methods
 A = GetProperty(H_hr,name)
end
methods (Static)
EQ_list = Test_TBSK_Var_gen(testmode)
end
methods(Hidden,Access= protected)
 H_hr = nn_sk_smart(H_hr,search_range,Accuracy,Rlength_cut)
 H_hr = nn_sk_sparse(H_hr,search_range,Accuracy,Rlength_cut)
 H_hr = H_TB_gen_SK(H_hr,options)
 H_hr = H_TB_gen_SK_sparse(H_hr,level_cut,para_filename,onsite_mode)
 Atom_smart_t = Atom_smart_t_gen(site1,site2)
 [nn_sparse_temp,Rnn_list] = nn_sparse_t_gen(site1,site2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
 [Rnn_list,nn_smart_t] = nn_smart_t_gen(Atom_smart_t,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
 hop = hop_gen(hop_pre,nn_level)
 [TBSK_hop,Coff] = TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
 [TBSK_hop,Coff] =  TBSK_hop_gen_sparse(site1,site2,Rij_cart,Rlength,nn_level)
end
end
