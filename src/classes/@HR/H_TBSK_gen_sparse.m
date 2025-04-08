function H_hr = H_TBSK_gen_sparse(H_hr,options)
arguments
H_hr HR;
options.level_cut {mustBeInteger} = -1;
options.onsite logical = false;
options.per_dir double = [1,1,1];
options.chiral logical = false;
options.spin logical = false;
options.method  = 'nn_sparse';
options.rough = false;
options.vectorL = [];
options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
options.para = struct();
options.Rd double = -1;
end
H_hr.num = true;
H_hr.coe = false;
if isempty(H_hr.nn_store)
if isempty(H_hr.nn_store_smart)
error('You have not run the function: nn/nn_smart');
else
if strcmp(options.method,'nn_sparse')
error('You should use the key value: ''method'',''nn_smart'' as input.');
end
end
end
N_orbit = H_hr.WAN_NUM;
if options.level_cut >0
select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
else
select_nn_store = H_hr.nn_store;
end
if options.chiral
i_L = select_nn_store(:,1);
j_L = select_nn_store(:,2);
element1_L = H_hr.elementL(i_L);
element2_L =  H_hr.elementL(j_L);
select_nn_store =select_nn_store(element1_L ~= element2_L,:);
end
if options.spin
i_L = select_nn_store(:,1);
j_L = select_nn_store(:,2);
spin1_L = H_hr.quantumL(i_L,4);
spin2_L = H_hr.quantumL(j_L,4);
select_nn_store =select_nn_store(spin1_L ~= spin2_L,:);
end
Rvector_L = select_nn_store(:,6:8);
[Rvector_L_unique,sorted_label,cut_slice] = TBkit.cut_tools(Rvector_L);
select_nn_store = select_nn_store(sorted_label,:);
i_L = select_nn_store(:,1);
j_L = select_nn_store(:,2);
Rlength_L = select_nn_store(:,9);
l_L = select_nn_store(:,3)./Rlength_L;
m_L = select_nn_store(:,4)./Rlength_L;
n_L = select_nn_store(:,5)./Rlength_L;
nn_level_L = select_nn_store(:,10);
L_1_L = H_hr.quantumL(i_L,2);
L_2_L = H_hr.quantumL(j_L,2);
m_1_L = H_hr.quantumL(i_L,3);
m_2_L = H_hr.quantumL(j_L,3);
NRPTS_ = size(cut_slice,1);
fprintf('Generating Numeric Hopping term ...')
if ~H_hr.overlap
NumHopping_L = HR.TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,...
'para',options.para,'Rd',options.Rd,...
'deltarule',options.deltarule ,...
'alpharule',options.alpharule ...
);
else
[NumHopping_L,NumOverlap_L] = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap',true);
end
fprintf('< ok >\n')
if strcmp(H_hr.Type,'mat')
H_hr = H_hr.add_empty_one(Rvector_L);
end
pb = TBkit_tool_outer.CmdLineProgressBar('Setting NRPT ');
for i = 1:NRPTS_
if ~H_hr.overlap
H_hr = H_hr.set_hop(...
NumHopping_L(cut_slice(i,1):cut_slice(i,2)),...
i_L(cut_slice(i,1):cut_slice(i,2)),...
j_L(cut_slice(i,1):cut_slice(i,2)),...
Rvector_L_unique(i,:),'set');
else
H_hr = H_hr.set_hop(...
NumHopping_L(cut_slice(i,1):cut_slice(i,2)),...
i_L(cut_slice(i,1):cut_slice(i,2)),...
j_L(cut_slice(i,1):cut_slice(i,2)),...
Rvector_L_unique(i,:),'set');
H_hr = H_hr.set_overlap(...
NumOverlap_L(cut_slice(i,1):cut_slice(i,2)),...
i_L(cut_slice(i,1):cut_slice(i,2)),...
j_L(cut_slice(i,1):cut_slice(i,2)),...
Rvector_L_unique(i,:),'set');
end
pb.print(i,NRPTS_,' ...');
end
if options.onsite
for i=1:N_orbit
onsite_sym_name = "E__"+string(H_hr.elementL(i,1))+"_"...
+string(H_hr.quantumL(i,2));...
H_hr = H_hr.set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'add');
end
end
if H_hr.overlap
for i=1:N_orbit
H_hr = H_hr.set_overlap_single(1,i,i,[0 0 0],'add');
end
end
pb.delete();
end
