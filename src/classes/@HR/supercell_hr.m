function H_hr = supercell_hr(H_hr,Ns,options)
arguments
H_hr HR;
Ns double = eye(H_hr.Dim);
options.Accuracy double = 1e-6;
options.force_list = false;
options.OBC = zeros(1,H_hr.Dim);
options.silence = false;
end
V= abs(round(det(Ns)));
OUT_WAN_NUM = H_hr.WAN_NUM*V ;
WANNUM= H_hr.WAN_NUM;
Accuracy = options.Accuracy;
pc_orb = H_hr.orbL;
[sc_orb,sc_vec,sc_elementL,sc_quantumL] = H_hr.supercell_orb(Ns,Accuracy);
OUT_H_hr = H_hr;
OUT_H_hr = OUT_H_hr.clean(OUT_WAN_NUM);
OUT_H_hr.orbL = sc_orb;
OUT_H_hr.quantumL = sc_quantumL;
OUT_H_hr.elementL = sc_elementL;
OUT_H_hr.Rm = Ns * OUT_H_hr.Rm;
NRPTS_  = H_hr.NRPTS;
if H_hr.overlap
NRPTS_S = size(H_hr.vectorL_overlap,1);
end
if ~options.silence
fprintf('Search done; begin to set hoppings\n');
fprintf('We may improve the perfomance later\n');
end
num_sc = size(sc_vec,1);
Accuracy_roundn = round(log(Accuracy)/log(10));
pc_orb_super = repmat(pc_orb,[num_sc,1]);
orb_sc_vevL = roundn(sc_orb*Ns - pc_orb_super,Accuracy_roundn);
Leak_L = ~ismembertol(orb_sc_vevL,sc_vec,Accuracy,'ByRows',true);
orb_sc_labelL = find(Leak_L,1);
if ~isempty(orb_sc_labelL) && strcmp(H_hr.Type,'mat') ||options.force_list
if ~options.silence
fprintf("supercell has leak sites, use list mode enforcely!");
end
H_hr = H_hr.rewrite();
OUT_H_hr = OUT_H_hr.rewrite();
end
switch H_hr.Type
case 'mat'
for icur_sc_vec = 1:num_sc
cur_sc_vec = double(sc_vec(icur_sc_vec,:));
if ~options.silence
pb = TBkit_tool_outer.CmdLineProgressBar(...
['Generate process: SUPERCELL(',...
num2str(icur_sc_vec),',',num2str(num_sc),') NRPT:']);
end
for ih = 1:NRPTS_
ind_R = double(H_hr.vectorL(ih,:));
indR_in_supercell=double(ind_R+cur_sc_vec)/Ns;
indR_in_supercell = roundn(indR_in_supercell,Accuracy_roundn);
sc_part=floor(indR_in_supercell);
orig_part=ind_R+cur_sc_vec-double(sc_part*Ns);
pair_ind = 9999;
for jcur_sc_vec = 1:num_sc
pair_sc_vec = sc_vec(jcur_sc_vec,:);
if pair_sc_vec==orig_part
if pair_ind ~= 9999
error("\n\nFound duplicate super cell vector!");
end
pair_ind=jcur_sc_vec;
end
end
if pair_ind==9999
disp(orig_part);
disp('Cant find sc in ');
disp(sc_vec);
disp('orig_part=ind_R+cur_sc_vec-sc_part*Ns;');
disp(ind_R);
disp(cur_sc_vec);
error("\n\nDid not find super cell vector!");
end
if H_hr.num
tmpHnum = H_hr.HnumL(:,:,ih);
Nonzero_list = find(tmpHnum ~= 0);
[hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
hi= hi + (icur_sc_vec-1)*H_hr.WAN_NUM ;
hj= hj + (pair_ind-1)*H_hr.WAN_NUM ;
tmp_mat = full(sparse(hi,hj,tmpHnum(Nonzero_list),OUT_WAN_NUM,OUT_WAN_NUM));
OUT_H_hr = OUT_H_hr.set_hop_mat(tmp_mat,sc_part,'add');
end
if H_hr.coe
tmpHcoe = H_hr.HcoeL(:,:,ih);
Nonzero_list = find(tmpHcoe ~= sym(0));
[hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
hi= hi + (icur_sc_vec-1)*H_hr.WAN_NUM ;
hj= hj + (pair_ind-1)*H_hr.WAN_NUM ;
tmp_mat = sym(zeros(OUT_WAN_NUM));
tmp_mat(sub2ind([OUT_WAN_NUM,OUT_WAN_NUM],hi,hj)) = tmpHcoe(Nonzero_list);
OUT_H_hr = OUT_H_hr.set_hop_mat(tmp_mat,sc_part,'symadd');
end
if ~options.silence
pb.print(ih,NRPTS_,' Hopping ...');
end
%                     fprintf("Generate process: SUPERCELL(%d,%d) NRPT(%d,%d) RUNINGTIME: %f s.\n",...
end
if ~options.silence
pb.delete();
end
if H_hr.overlap
pb = TBkit_tool_outer.CmdLineProgressBar(...
['Generate process: SUPERCELL(S mat)(',...
num2str(icur_sc_vec),',',num2str(num_sc),') NRPT:']);
for is = 1:NRPTS_S
ind_R_S = double(H_hr.vectorL_overlap(is,:));
indR_in_supercell=double(ind_R_S+cur_sc_vec)/Ns;
indR_in_supercell = roundn(indR_in_supercell,Accuracy_roundn);
sc_part=floor(indR_in_supercell);
orig_part=ind_R_S+cur_sc_vec-double(sc_part*Ns);
pair_ind = 9999;
for jcur_sc_vec = 1:num_sc
pair_sc_vec = sc_vec(jcur_sc_vec,:);
if pair_sc_vec==orig_part
if pair_ind ~= 9999
error("\n\nFound duplicate super cell vector!");
end
pair_ind=jcur_sc_vec;
end
end
if pair_ind==9999
disp(orig_part);
disp('Cant find sc in ');
disp(sc_vec);
disp('orig_part=ind_R+cur_sc_vec-sc_part*Ns;');
disp(ind_R);
disp(cur_sc_vec);
error("\n\nDid not find super cell vector!");
end
if H_hr.num
tmpSnum = H_hr.SnumL(:,:,is);
Nonzero_list = find(tmpSnum ~= 0);
[si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
si= si + (icur_sc_vec-1)*H_hr.WAN_NUM ;
sj= sj + (pair_ind-1)*H_hr.WAN_NUM ;
tmp_mat = full(sparse(si,sj,tmpSnum(Nonzero_list),OUT_WAN_NUM,OUT_WAN_NUM));
OUT_H_hr = OUT_H_hr.set_overlap_mat(tmp_mat,sc_part,'add');
end
if H_hr.coe
tmpScoe = H_hr.ScoeL(:,:,is);
Nonzero_list = find(tmpScoe ~= sym(0));
[si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
si= si + (icur_sc_vec-1)*H_hr.WAN_NUM ;
sj= sj + (pair_ind-1)*H_hr.WAN_NUM ;
tmp_mat = sym(zeros(OUT_WAN_NUM));
tmp_mat(sub2ind([OUT_WAN_NUM,OUT_WAN_NUM],si,sj)) = tmpScoe(Nonzero_list);
OUT_H_hr = OUT_H_hr.set_overlap_mat(tmp_mat,sc_part,'symadd');
end
pb.print(is,NRPTS_S,' Overlap ...');
end
pb.delete();
%                     fprintf("Generate process: SUPERCELL(%d,%d) NRPT(%d,%d) RUNINGTIME: %f s.\n",...
end
end
case 'list'
if H_hr.num
HnumList = repmat(H_hr.HnumL,[num_sc,1]);
nHopping = length(H_hr.HnumL);
end
if H_hr.coe
HcoeList = repmat(H_hr.HcoeL,[num_sc,1]);
nHopping = length(H_hr.HcoeL);
end
VectorList = double(H_hr.vectorL);
OutVectorList = repmat(VectorList,[num_sc,1]);
if ~options.silence
pb = TBkit_tool_outer.CmdLineProgressBar(...
'Generate process: SUPERCELL:');
end
for icur_sc_vec = 1:num_sc
cur_sc_vec = double(sc_vec(icur_sc_vec,:));
hiL = VectorList(:,H_hr.Dim+1);
hjL = VectorList(:,H_hr.Dim+2);
ind_RL = double(VectorList(:,1:H_hr.Dim));
ind_R_in_supercellL = double(ind_RL+cur_sc_vec)/Ns;
ind_R_in_supercellL = roundn(ind_R_in_supercellL,Accuracy_roundn);
sc_partL=floor(ind_R_in_supercellL);
orig_partL=ind_RL+cur_sc_vec-double(sc_partL*Ns);
[~,pair_indL] = ismember(orig_partL,sc_vec,'rows');
sc_hjL = hjL+(pair_indL-1)*WANNUM;
sc_hiL = hiL + (icur_sc_vec-1)*WANNUM;
indRtiL = pc_orb(hiL,:);
indRti_in_supercellL = sc_orb(sc_hiL,:);
real_sc_vecL = indRti_in_supercellL*Ns - indRtiL;
real_sc_vecL = round(real_sc_vecL);
indRtjL = real_sc_vecL+ind_RL+pc_orb(hjL,:);
indRtj_in_supercellL = roundn(double(indRtjL)/Ns,Accuracy_roundn);
indR_in_supercellL  = floor(indRtj_in_supercellL);
OutVectorList((icur_sc_vec-1)*nHopping+1:(icur_sc_vec)*nHopping,:) = ...
[indR_in_supercellL,sc_hiL,sc_hjL];
if ~options.silence
pb.print(icur_sc_vec,num_sc,' ...');
end
end
if ~options.silence
pb.delete();
end
if H_hr.num
OUT_H_hr.HnumL = HnumList;
end
if H_hr.coe
OUT_H_hr.HcoeL = HcoeList;
end
OUT_H_hr.vectorL = OutVectorList;
end
H_hr = OUT_H_hr;
H_hr.Basis_num = OUT_WAN_NUM;
H_hr = OpenBoundary(H_hr,options.OBC);
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
