function H_hr = init(H_hr,options)
arguments
H_hr HR;
options.level_cut {mustBeInteger} = 1;
options.level_list = [];
options.onsite logical = false;
options.per_dir double = [1,1,1];
options.chiral logical = false;
options.spin logical = false;
options.method  = 'nn_sparse';
options.rough logical= false;
options.hermitize = true;
options.vectorL = [];
options.fast logical= false;
end
if ~isempty(options.vectorL)
switch size(options.vectorL,2)
case 3
case 5
end
if options.hermitize
H_hr = H_hr.hermitize();
end
return;
end
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
N_tot = H_hr.WAN_NUM;
level_cut = options.level_cut;
if ~strcmp(H_hr.Type,'list')
H_hr = H_hr.rewrite();
end
H_hr.coe = true;
H_hr.num = false;
if isempty(options.level_list)
select_nn_store = H_hr.nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
else
select_nn_store = H_hr.nn_store(ismember(H_hr.nn_store(:,10),options.level_list),:);
end
if options.fast
H_hr.vectorhopping = true;
if options.chiral
element1L =  H_hr.elementL(select_nn_store(:,1));
element2L =  H_hr.elementL(select_nn_store(:,2));
select_nn_store(element1L == element2L,:) = [];
end
if options.spin
spin1L = H_hr.quantumL(select_nn_store(:,1),4);
spin2L = H_hr.quantumL(select_nn_store(:,2),4);
select_nn_store(spin1L == spin2L,:) = [];
end
if options.onsite
onsite_nnL = [(1:N_orbit)',(1:N_orbit)',zeros(N_orbit,8)];
select_nn_store = [select_nn_store;onsite_nnL];
end
nselect_nn_store = size(select_nn_store,1);
H_hr.vectorL = select_nn_store(:,[6,7,8,1,2]);
H_hr.AvectorL = eye(nselect_nn_store);
H_hr.BvectorL = eye(nselect_nn_store);
H_hr.CvectorL = eye(nselect_nn_store*2);
else
if options.onsite
if strcmp(options.method,'nn_smart')
for i=1:N_orbit
A = sym(['A','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
B = sym(['B','_',num2str(i),'_',num2str(i),'_0_ubar'],'real');
H_hr = H_hr.set_hop_single( A+1i*B,i,i,[0 0 0],'sym');
end
else
onsite_nnL = [(1:N_orbit)',(1:N_orbit)',zeros(N_orbit,8)];
H_hr.nn_store = [H_hr.nn_store;onsite_nnL];
end
end
if strcmp(options.method,'nn_smart')
pb = TBkit_tool_outer.CmdLineProgressBar('Setting ');
for i=1:N_orbit
%fprintf("setting (%4d/%4d) th orbital ... \n",i,N_orbit);
for j=1:N_tot
set_or_not = true;
if options.chiral
element1 =  H_hr.elementL(i);
element2 =  H_hr.elementL(j);
if element1 == element2
set_or_not = false;
end
end
if options.spin
spin1 = H_hr.quantumL(i,4);
spin2 = H_hr.quantumL(j,4);
if spin1 == spin2
set_or_not = false;
end
end
if set_or_not
nn = H_hr.nn_store_smart(i,j).nn;
for k = 1:size(nn,1)
if nn(k).nn_level <= level_cut
vector = nn(k).R_vector.*per_dir;
A = HR.SymbolicVarible("A",vector,[i,j],nn(k).nn_level);
B = HR.SymbolicVarible("B",vector,[i,j],nn(k).nn_level);
SymHopping = A+1i*B;
H_hr = H_hr.set_hop(SymHopping,i,j,vector,'sym');
end
end
end
pb.print([i,j],[N_orbit,N_tot],' th orbital ...');
end
end
pb.delete();
else
pb = TBkit_tool_outer.CmdLineProgressBar('Setting ');
if options.chiral
element1L =  H_hr.elementL(select_nn_store(:,1));
element2L =  H_hr.elementL(select_nn_store(:,2));
select_nn_store(element1L == element2L,:) = [];
end
if options.spin
spin1L = H_hr.quantumL(select_nn_store(:,1),4);
spin2L = H_hr.quantumL(select_nn_store(:,1),4);
select_nn_store(spin1L == spin2L,:) = [];
end
nselect_nn_store = size(select_nn_store,1);
iL = select_nn_store(:,1);STRiL = string(iL);
jL = select_nn_store(:,2);STRjL = string(jL);
R1L = select_nn_store(:,6);minusR1L = R1L <0;STRR1L = string(abs(R1L));STRR1L(minusR1L) =STRR1L(minusR1L)+"_bar";
R2L = select_nn_store(:,7);minusR2L = R2L <0;STRR2L = string(abs(R2L));STRR2L(minusR2L) =STRR2L(minusR2L)+"_bar";
R3L = select_nn_store(:,8);minusR3L = R3L <0;STRR3L = string(abs(R3L));STRR3L(minusR3L) =STRR3L(minusR3L)+"_bar";
nn_levelL = select_nn_store(:,10);STRnnL = string(nn_levelL)+"_ubar";
SuperscriptL = repmat("__",[nselect_nn_store 1]);
SubscriptL = repmat("_",[nselect_nn_store 1]);
affix_L = SuperscriptL+STRR1L+SuperscriptL+STRR2L+SuperscriptL+STRR3L+...
SubscriptL+STRiL+SubscriptL+STRjL+SubscriptL+STRnnL;
AL = "A" + affix_L;
BL = "B" + affix_L;
AsymL = str2sym(AL);
BsymL = str2sym(BL);
assume(AsymL,'real');
assume(BsymL,'real');
H_hr.vectorL = select_nn_store(:,[6,7,8,1,2]);
H_hr.HcoeL = AsymL+1i*BsymL;
pb.delete();
end
end
if options.hermitize
H_hr = H_hr.hermitize();
end
end
