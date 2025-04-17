function H_hr = unfold_hr(H_hr,Ns,options)
%UNFOLD_HR Unfold supercell Hamiltonian to primitive cell
%
%   H_HR = UNFOLD_HR(H_HR, NS, OPTIONS) transforms supercell Hamiltonian
%   to primitive cell representation
%
%   Inputs:
%       H_hr       - HR object in supercell representation
%       Ns         - Supercell transformation matrix (3×3)
%       options.Accuracy    - Numerical tolerance (default = 1e-6)
%       options.force_list  - Force list storage format (default = false)
%       options.orb_idL     - Optional orbital ID list
%
%   Output:
%       H_hr - HR object in primitive cell representation
%
%   Note:
%       Requires det(Ns) to be integer and positive
%
%   See also HR.FOLD_HR, UNFOLD_ORB
arguments
H_hr HR;
Ns double = eye(3);
options.Accuracy double = 1e-6;
options.force_list = false;
options.orb_idL = [];
end
V= 1/abs(round(det(Ns)));
OUT_WAN_NUM = H_hr.WAN_NUM*V;
if  rem(OUT_WAN_NUM,1)
error('The supercell matrix for primitive cell is not right.\n');
else
end
WANNUM= H_hr.WAN_NUM;
Accuracy = options.Accuracy;
sc_orbL = H_hr.orbL;
[pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,sc_orb_idL,~,pc_orb_selectL] = H_hr.unfold_orb(Ns,Accuracy,options.orb_idL);
OUT_H_hr = H_hr;
OUT_H_hr = OUT_H_hr.clean(WANNUM);
OUT_H_hr.orbL = pc_orb;
OUT_H_hr.quantumL = pc_quantumL;
OUT_H_hr.elementL = pc_elementL;
OUT_H_hr.Rm = Ns\OUT_H_hr.Rm;
NRPTS_  = H_hr.NRPTS;
if strcmp(H_hr.Type,'sparse')
H_hr = H_hr.full();
end
fprintf('Search done; begin set hoppings\n');
fprintf('We can improve the perfomance later\n');
Accuracy_roundn = round(log(Accuracy)/log(10));
if (strcmp(H_hr.Type,'mat')) ||options.force_list
fprintf("Attention enforce list mode in the folding process !\n");
H_hr = H_hr.rewrite();
OUT_H_hr = OUT_H_hr.rewrite();
end
switch H_hr.Type
case 'list'
VectorList = double(H_hr.vectorL);
pb = CmdLineProgressBar(...
'Generate process: UNFOLDING:\n');
sc_hiL = VectorList(:,H_hr.Dim+1);
sc_hjL = VectorList(:,H_hr.Dim+2);
hjL_orb_id_in_primitiveL = sc_orb_idL(sc_hjL);
hiL_orb_id_in_primitiveL = sc_orb_idL(sc_hiL);
Npc_orb_selectL = find(pc_orb_selectL);
SelectedL = ismember(sc_hiL,Npc_orb_selectL);
if H_hr.num
HnumList = H_hr.HnumL(SelectedL,:);
end
if H_hr.coe
HcoeList = H_hr.HcoeL(SelectedL,:);
end
hjL_orb_id_in_primitiveL = hjL_orb_id_in_primitiveL(SelectedL);
hiL_orb_id_in_primitiveL = hiL_orb_id_in_primitiveL(SelectedL);
Selected_vectorL = VectorList(SelectedL,:);
ind_R_in_supercellL = double(Selected_vectorL(:,1:H_hr.Dim));
Selected_sc_hiL = Selected_vectorL(:,H_hr.Dim+1);
Selected_sc_hjL = Selected_vectorL(:,H_hr.Dim+2);
TijL_in_supercellL = ind_R_in_supercellL + ...
sc_orbL(Selected_sc_hjL,:) - sc_orbL(Selected_sc_hiL,:);
TijL_in_primitiveL = TijL_in_supercellL*Ns;
hiL_orbL_in_primitiveL = pc_orbL_full(Selected_sc_hiL,:);
hjL_plusR_orbL_in_primitiveL = hiL_orbL_in_primitiveL + TijL_in_primitiveL;
RvectorL_in_primitiveL = floor(hjL_plusR_orbL_in_primitiveL);
OutVectorList = [RvectorL_in_primitiveL,hiL_orb_id_in_primitiveL.',hjL_orb_id_in_primitiveL.'];
pb.delete();
if H_hr.num
OUT_H_hr.HnumL = HnumList;
end
if H_hr.coe
OUT_H_hr.HcoeL = HcoeList;
end
OUT_H_hr.vectorL = OutVectorList;
otherwise
end
H_hr = OUT_H_hr;
H_hr.Basis_num = OUT_WAN_NUM;
if options.force_list
if strcmp(H_hr.Type,'mat')
H_hr = H_hr.rewrite();
end
else
if ~strcmp(H_hr.Type,'mat')
H_hr = H_hr.rewind();
end
end
end
