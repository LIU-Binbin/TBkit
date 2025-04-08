function H_hr = cut_piece(H_hr,repeatnum,fin_dir,glue_edges,vacuum_mode)
arguments
H_hr HR;
repeatnum double{mustBeInteger} =  10;
fin_dir double{mustBeMember(fin_dir,[1,2,3,4])} = 3;
glue_edges logical = false;
vacuum_mode logical = false;
end
Ns = [1 0 0;0 1 0;0 0 1];
Ns(fin_dir,:) = Ns(fin_dir,:) * repeatnum;
fin_dir_list = [0 0 0];
[sc_orb,~,sc_elementL,sc_quantumL] = H_hr.supercell_orb(Ns);
sc_orb = double(sc_orb);
if vacuum_mode
fin_dir_list(fin_dir) = 1;
end
if isempty(H_hr.sites)
if exist('POSCAR','file')
[Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,~]=HR.POSCAR_read('POSCAR','vasp');
H_hr = H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
else
end
else
try
H_hr = H_hr.supercell(Ns,'POSCAR_super_fin',H_hr.Rm,H_hr.sites,H_hr.Atom_name,H_hr.Atom_num,fin_dir_list);
catch
end
end
norb =size(H_hr.orbL,1);
if norb ~= H_hr.WAN_NUM
error("\n\nOribital_init is wrong,please give a right orbital init or just forget this parm!");
end
if repeatnum<1
error("\n\nArgument num must be positive!");
end
if repeatnum == 1 && glue_edges
error("\n\nCan't have num==1 and glueing of the edges!");
end
if vacuum_mode
Rmlength1 = norm (H_hr.Rm(1,:));
Rmlength2 = norm (H_hr.Rm(2,:));
Rmlength3 = norm (H_hr.Rm(3,:));
Rm_s_fin_add = [10*H_hr.Rm(1,:)*fin_dir_list(1)/Rmlength1;...
10*H_hr.Rm(2,:)*fin_dir_list(2)/Rmlength2;...
10*H_hr.Rm(3,:)*fin_dir_list(3)/Rmlength3];
Rm_s_fin = H_hr.Rm + Rm_s_fin_add ;
Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
[nfinorb,~ ]= size(sc_orb);
for  i = 1:nfinorb
Rr_orb = sc_orb(i,:)*H_hr.Rm;
Rr_s_fin = Rr_orb + Rr_s_fin_add;
Rc_s_fin = Rr_s_fin /Rm_s_fin;
sc_orb(i,:) = Rc_s_fin ;
end
end
OUT_WAN_NUM = H_hr.WAN_NUM*repeatnum ;
vectorList = double(H_hr.vectorL);
NRPT_seq = false(H_hr.NRPTS,1);
if H_hr.overlap
vectorList_overlap = double(H_hr.vectorL_overlap);
NRPTS_S = size(vectorList_overlap,1);
NRPT_seq_S = false(NRPTS_S,1);
end
WANNUM = H_hr.WAN_NUM;
if strcmp(H_hr.Type,'sparse')
OUT_HnumL{H_hr.NRPTS} = sparse(OUT_WAN_NUM,OUT_WAN_NUM);
for i =1: H_hr.NRPTS-1
OUT_HnumL{i} = sparse(OUT_WAN_NUM,OUT_WAN_NUM);
end
pb = TBkit_tool_outer.CmdLineProgressBar('NRPT : ');
for iN = 1:H_hr.NRPTS
ind_R = vectorList(iN,:);
jump_fin=ind_R(fin_dir);
pb.print(iN,H_hr.NRPTS);
[ilist,jlist,amplist] = find(H_hr.HnumL{iN});
nhopping = length(amplist);
for ih = 1:nhopping
i = ilist(ih);
j = jlist(ih);
amp = amplist(ih);
for icur_sc_vec = 1:repeatnum
hi= i + (icur_sc_vec-1)*WANNUM ;
hj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
to_add=1;
if ~glue_edges
if hj <= 0 || hj > OUT_WAN_NUM
to_add=0;
end
else
hj= mod(hj,OUT_WAN_NUM);
if hj ==0
hj = OUT_WAN_NUM;
end
end
if to_add == 1
OUT_HnumL{iN}(hi,hj) = amp;
end
end
end
end
pb.delete();
H_hr.HnumL = OUT_HnumL;
H_hr.orbL = sc_orb    ;
H_hr.Type = 'sparse'     ;
elseif strcmp(H_hr.Type,'mat')
OUT_HnumL = zeros(OUT_WAN_NUM,OUT_WAN_NUM,H_hr.NRPTS);
NRPTS_record = H_hr.NRPTS;
if H_hr.overlap
OUT_SnumL = zeros(OUT_WAN_NUM,OUT_WAN_NUM,H_hr.NRPTS);
end
pb = TBkit_tool_outer.CmdLineProgressBar('NRPT : ');
for ih = 1:H_hr.NRPTS
ind_R = vectorList(ih,:);
jump_fin=ind_R(fin_dir);
tmpHnum = H_hr.HnumL(:,:,ih);
Nonzero_list = find(tmpHnum ~= 0);
[hi,hj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
pb.print(ih,H_hr.NRPTS);
HI_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
HJ_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
CONTAIN_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
rm_count  = 1;
for icur_sc_vec = 1:repeatnum
ind_R_tmp = ind_R;
hi_list = hi + (icur_sc_vec-1)*WANNUM ;
hj_list = hj + (icur_sc_vec+jump_fin-1)*WANNUM ;
if ~glue_edges
Contain_list = find(hj_list > 0 & hj_list <= OUT_WAN_NUM);
hj_list = hj_list(Contain_list);
hi_list = hi_list(Contain_list);
ind_R_tmp(fin_dir) = 0;
else
hj_list= mod(hj_list,OUT_WAN_NUM);
ind_R_tmp(fin_dir) = floor(hj./OUT_WAN_NUM);
hj_list(hj_list ==0) = OUT_WAN_NUM;
error('There is a bug, contact me: parkman@buaa.edu.cn');
end
[~,IH] = ismember(ind_R_tmp,vectorList,'rows');
nContain_list = length(Contain_list) ;
HI_list(rm_count:rm_count +nContain_list-1) = hi_list ;
HJ_list(rm_count:rm_count +nContain_list-1 ) = hj_list;
CONTAIN_list(rm_count:rm_count +nContain_list-1 ) = Contain_list;
rm_count = rm_count + nContain_list;
end
HI_list(rm_count:end) = [];
HJ_list(rm_count:end) = [];
CONTAIN_list(rm_count:end) = [];
if IH == 0
IH = NRPTS_record+1;
NRPTS_record = NRPTS_record+1;
vectorList(IH,:) = ind_R_tmp;
OUT_HnumL(:,:,IH) = full(sparse(HI_list,HJ_list,tmpHnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
else
OUT_HnumL(:,:,IH) = OUT_HnumL(:,:,IH) + full(sparse(HI_list,HJ_list,tmpHnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
end
NRPT_seq(IH)  = true;
end
pb.delete();
if H_hr.overlap
pb = TBkit_tool_outer.CmdLineProgressBar('NRPT_S : ');
for is = 1:NRPTS_S
ind_R_S = vectorList_overlap(is,:);
jump_fin=ind_R_S(fin_dir);
tmpSnum = H_hr.SnumL(:,:,is);
Nonzero_list = find(tmpSnum ~= 0);
[si,sj] = ind2sub([WANNUM,WANNUM],Nonzero_list);
pb.print(is,NRPTS_S);
SI_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
SJ_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
CONTAIN_list = zeros(1,OUT_WAN_NUM*OUT_WAN_NUM);
rm_count  = 1;
for icur_sc_vec = 1:repeatnum
ind_R_tmp  = ind_R_S;
si_list = si + (icur_sc_vec-1)*WANNUM ;
sj_list = sj + (icur_sc_vec+jump_fin-1)*WANNUM ;
if ~glue_edges
Contain_list = find(sj_list > 0 & sj_list <= OUT_WAN_NUM);
sj_list = sj_list(Contain_list) ;
si_list = si_list(Contain_list) ;
ind_R_tmp(fin_dir) = 0;
else
sj_list= mod(sj_list,OUT_WAN_NUM);
ind_R_tmp(fin_dir) = floor(hj./OUT_WAN_NUM);
sj_list(sj_list ==0) = OUT_WAN_NUM;
error('There is a bug, contact me: parkman@buaa.edu.cn');
end
[~,IS] = ismember(ind_R_tmp,vectorList_overlap,'rows');
NRPT_seq_S(IS) = true;
nContain_list = length(Contain_list) ;
SI_list(rm_count:rm_count +nContain_list-1) = si_list ;
SJ_list(rm_count:rm_count +nContain_list-1 ) = sj_list;
CONTAIN_list(rm_count:rm_count +nContain_list-1 ) = Contain_list;
rm_count = rm_count + nContain_list;
end
SI_list(rm_count:end) = [];
SJ_list(rm_count:end) = [];
CONTAIN_list(rm_count:end) = [];
OUT_SnumL(:,:,IS)= OUT_SnumL(:,:,IS)+full(sparse(SI_list,SJ_list,tmpSnum(Nonzero_list(CONTAIN_list)),OUT_WAN_NUM,OUT_WAN_NUM));
end
end
H_hr.vectorL = int8(vectorList);
H_hr.HnumL = OUT_HnumL;
if H_hr.overlap
H_hr.SnumL = OUT_SnumL;
end
if H_hr.overlap
H_hr = H_hr.reseq(':',NRPT_seq,NRPT_seq_S);
else
H_hr = H_hr.reseq(':',NRPT_seq);
end
elseif strcmp(H_hr.Type,'list')
OUT_HnumL = zeros(H_hr.NRPTS*repeatnum,1);
vector_list = H_hr.vectorL;
OUT_vectorList = zeros(H_hr.NRPTS*repeatnum,H_hr.Dim+2);
pb = TBkit_tool_outer.CmdLineProgressBar('H:NRPT : ');
count = 0;
for ih = 1:H_hr.NRPTS
ind_R = vectorList(ih,1:H_hr.Dim);
jump_fin=ind_R(fin_dir);
pb.print(ih,H_hr.NRPTS);
i = vector_list(ih,H_hr.Dim+1);
j = vector_list(ih,H_hr.Dim+2);
amp = H_hr.HnumL(ih);
if norm(amp) > 0
for icur_sc_vec = 1:repeatnum
hi= i + (icur_sc_vec-1)*WANNUM ;
hj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
to_add=1;
if ~glue_edges
if hj <= 0 || hj > OUT_WAN_NUM
to_add=0;
end
else
hj= mod(hj,OUT_WAN_NUM);
if hj ==0
hj = OUT_WAN_NUM;
end
end
if to_add == 1
count = count+1;
OUT_HnumL(count) = amp;
OUT_vectorList(count,:) = [ind_R,hi,hj];
end
end
end
end
pb.delete();
fprintf('remove %d empty orbs',H_hr.NRPTS*repeatnum- count);
OUT_HnumL(count+1:end,:) = [];
OUT_vectorList(count+1:end,:) = [];
H_hr.HnumL = OUT_HnumL;
OUT_vectorList(:,fin_dir) = 0;
H_hr.vectorL = OUT_vectorList;
if H_hr.overlap
NRPTS_S = size(H_hr.vectorL_overlap,2);
vectorList_overlap = H_hr.vectorL_overlap;
OUT_vectorList_overlap = zeros(NRPTS_S*repeatnum,H_hr.Dim+2);
OUT_SnumL = zeros(NRPTS_S*repeatnum,1);
pb = TBkit_tool_outer.CmdLineProgressBar('S:NRPT : ');
count = 0;
for ih = 1:NRPTS_S
ind_R = vectorList_overlap(ih,1:H_hr.Dim);
jump_fin=ind_R(fin_dir);
pb.print(ih,H_hr.NRPTS);
i = vectorList_overlap(ih,H_hr.Dim+1);
j = vectorList_overlap(ih,H_hr.Dim+2);
amp = H_hr.SnumL(ih);
if norm(amp) > 0
for icur_sc_vec = 1:repeatnum
si= i + (icur_sc_vec-1)*WANNUM ;
sj= j + (icur_sc_vec+jump_fin-1)*WANNUM ;
to_add=1;
if ~glue_edges
if sj <= 0 || sj > OUT_WAN_NUM
to_add=0;
end
else
sj= mod(sj,OUT_WAN_NUM);
if sj ==0
sj = OUT_WAN_NUM;
end
end
if to_add == 1
count = count+1;
OUT_SnumL(count) = amp;
OUT_vectorList_overlap(count,:) = [ind_R,si,sj];
end
end
end
end
pb.delete();
H_hr.SnumL = OUT_SnumL;
OUT_vectorList_overlap(:,fin_dir) = 0;
H_hr.vectorL_overlap = OUT_vectorList_overlap;
end
end
H_hr.num = true    ;
H_hr.coe = false    ;
H_hr.orbL = sc_orb  ;
H_hr.quantumL = sc_quantumL;
H_hr.elementL = sc_elementL;
end
