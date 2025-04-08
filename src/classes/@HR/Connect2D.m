function H_hr_n = Connect2D(H_hr_n,H_hr_unitcell,opt)
arguments
H_hr_n HR;
H_hr_unitcell HR;
opt.Sequ double = [1:20];
opt.findir double = 2;
opt.Subsequ {mustBeMember(opt.Subsequ,{'normal','reverse','combined'})} = 'normal';
opt.reverse_sequ double = [H_hr_unitcell.WAN_NUM:-1:1];
opt.combined_norm double = [2:19];
end
WANNUM = H_hr_n.WAN_NUM;
WANNUM_unit = H_hr_unitcell.WAN_NUM;
if opt.findir==2
Nrep_y = length(opt.Sequ);
Nrep_x = WANNUM/WANNUM_unit/Nrep_y;
Nrep_tot = Nrep_x*Nrep_y;
for k = -1:1:1
slct = ismember(H_hr_unitcell.vectorL(:,1:H_hr.Dim),[1 k 0],'rows');
if all(slct==0)
continue;
end
if strcmp(opt.Subsequ,'normal')
[row,col,num] = find(H_hr_unitcell.HnumL(:,:,slct) );
elseif strcmp(opt.Subsequ,'reverse')
HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ;
HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
[row,col,num] = find(HnumL_reverse );
disp('This reverse mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
if num~=flip(num)
warning('This unitcell HR may not be suitable for the reverse mode. Check the relection symmetry for Hr_unitcell!');
end
elseif strcmp(opt.Subsequ,'combined')
HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ;
[row_norm,col_norm,num_norm] = find(HnumL_tmp );
HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
[row_rev,col_rev,num_rev] = find(HnumL_reverse );
disp('This combined mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
if num_rev~=flip(num_rev)
warning('This unitcell HR may not be suitable for the combined mode. Check the relection symmetry for Hr_unitcell!');
end
if max(opt.combined_norm)>WANNUM || min(opt.combined_norm)<1
warning('Check the combined_norm range!');
end
end
if strcmp(opt.Subsequ,'combined')
count = 1;norm_length = length(opt.combined_norm);
for i = 1:length(opt.Sequ)
labR = Nrep_y*(Nrep_x-1)+i;
labL =  opt.Sequ(i) +k;
if labL>0&&labL<=Nrep_y
if count<=norm_length && i == opt.combined_norm(count)
for j =1: length(row_norm)
H_hr_n = H_hr_n.set_hop(...
num_norm(j),WANNUM_unit*(labR-1)+row_norm(j),WANNUM_unit*(labL-1)+col_norm(j),[1 0 0],'add');
H_hr_n = H_hr_n.set_hop(...
conj(num_norm(j)),WANNUM_unit*(labL-1)+col_norm(j),WANNUM_unit*(labR-1)+row_norm(j),[-1 0 0],'add');
end
count = count +1;
else
for j =1: length(row_rev)
H_hr_n = H_hr_n.set_hop(...
num_rev(j),WANNUM_unit*(labR-1)+row_rev(j),WANNUM_unit*(labL-1)+col_rev(j),[1 0 0],'add');
H_hr_n = H_hr_n.set_hop(...
conj(num_rev(j)),WANNUM_unit*(labL-1)+col_rev(j),WANNUM_unit*(labR-1)+row_rev(j),[-1 0 0],'add');
end
end
end
end
else
for i = 1:length(opt.Sequ)
labR = Nrep_y*(Nrep_x-1)+i;
labL =  opt.Sequ(i) +k;
if labL>0&&labL<=Nrep_y
for j =1: length(row)
H_hr_n = H_hr_n.set_hop(...
num(j),WANNUM_unit*(labR-1)+row(j),WANNUM_unit*(labL-1)+col(j),[1 0 0],'add');
H_hr_n = H_hr_n.set_hop(...
conj(num(j)),WANNUM_unit*(labL-1)+col(j),WANNUM_unit*(labR-1)+row(j),[-1 0 0],'add');
end
end
end
end
end
elseif opt.findir==1
Nrep_x = length(opt.Sequ);
Nrep_y = WANNUM/WANNUM_unit/Nrep_x;
for k = -1:1:1
slct = ismember(H_hr_unitcell.vectorL(:,1:H_hr.Dim),[k 1 0],'rows');
if all(slct==0)
continue;
end
if strcmp(opt.Subsequ,'normal')
[row,col,num] = find(H_hr_unitcell.HnumL(:,:,slct) );
elseif strcmp(opt.Subsequ,'reverse')
HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct) ;
HnumL_reverse = HnumL_tmp(:,opt.reverse_sequ);
[row,col,num] = find(HnumL_reverse );
disp('This reverse mode only work for reflection symmetric inter-cell hoppings, and may have bugs for double layer materials');
if num~=flip(num)
warning('This unitcell HR may not be suitable for the reverse mode. Check the relection symmetry for Hr_unitcell!');
end
end
for i = 1:length(opt.Sequ)
labR = i*(Nrep_y);
tmp = (opt.Sequ(i)-1+k);
labL = Nrep_y*tmp + 1;
if tmp>=0&&tmp<Nrep_x
for j =1: length(row)
H_hr_n = H_hr_n.set_hop(...
num(j),WANNUM_unit*(labR-1)+row(j),WANNUM_unit*(labL-1)+col(j),[0 1 0],'add');
H_hr_n = H_hr_n.set_hop(...
conj(num(j)),WANNUM_unit*(labL-1)+col(j),WANNUM_unit*(labR-1)+row(j),[0 -1 0],'add');
end
end
end
end
end
end
