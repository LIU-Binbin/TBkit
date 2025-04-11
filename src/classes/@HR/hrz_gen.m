function H_hr= hrz_gen(H_hr,kpoints_f,fin_dir,mode)
%HRZ_GEN Generate horizontal representation of HR object
%   Transforms an HR object into a horizontal representation based on
%   specified k-points and direction.
%
%   Syntax:
%       H_hr = hrz_gen(H_hr)
%       H_hr = hrz_gen(H_hr,kpoints_f)
%       H_hr = hrz_gen(H_hr,kpoints_f,fin_dir)
%       H_hr = hrz_gen(H_hr,kpoints_f,fin_dir,mode)
%
%   Inputs:
%       H_hr - Input HR object
%       kpoints_f - k-points in fractional coordinates (default: [0,0,0])
%       fin_dir - Direction for horizontalization (default: 3)
%       mode - Mode of operation ('0D', '1D', or '2D') (default: '1D')
%
%   Outputs:
%       H_hr - Modified HR object in horizontal representation
%
%   Example:
%       H_hr_1D = hrz_gen(H_hr, [0 0 0], 1, '1D');

if nargin < 4
mode = '1D';
end
if nargin <3
fin_dir = 3;
end
if nargin < 2
kpoints_f = [0 ,0 ,0];
end
import TBkit_tool.*
factor_list = exp(1i*2*pi*(double(H_hr.vectorL)*kpoints_f.'));
[Num,Coe] = H_hr.NumOrCoe;
switch H_hr.Type
case {'sparse'}
Hnum_kf_list = reshape(full(cell2mat(H_hr.HnumL)),H_hr.WAN_NUM,H_hr.WAN_NUM,H_hr.NRPTS);
case 'mat'
if Num
Hnum_kf_list = H_hr.HnumL;
end
if Coe
Hcoe_kf_list = H_hr.HcoeL;
end
case 'list'
otherwise
end
switch H_hr.Type
case {'sparse'}
for i = 1 : H_hr.NRPTS
Hnum_kf_list(:,:,i) = Hnum_kf_list(:,:,i)*factor_list(i);
end
case 'mat'
if Num
for i = 1 : H_hr.NRPTS
Hnum_kf_list(:,:,i) = Hnum_kf_list(:,:,i)*factor_list(i);
end
end
if Coe
for i = 1 : H_hr.NRPTS
Hcoe_kf_list(:,:,i) = Hcoe_kf_list(:,:,i)*factor_list(i);
end
end
end
switch mode
case '0D'
switch H_hr.Type
case {'sparse'}
H_hr.HnumL = sum(Hnum_kf_list,3);
case {'mat'}
if Num
H_hr.HnumL = sum(Hnum_kf_list,3);
end
if Coe
H_hr.HcoeL = sum(Hcoe_kf_list,3);
end
end
case '1D'
[vector_list_new,sort_label] = sortrows(H_hr.vectorL,fin_dir) ;
switch H_hr.Type
case {'sparse'}
Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label);
case {'mat'}
if Num
Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label);
end
if Coe
Hcoe_kf_list_new = Hcoe_kf_list(:,:,sort_label);
end
end
[unique_dir,unique_label]= unique(vector_list_new(:,fin_dir));
cutlist(:,1)= unique_label;
cutlist(1:end-1,2)= unique_label(2:end)-1;
cutlist(end,2) = H_hr.NRPTS;
vector_list_hrz = zeros(length(unique_dir),3);
NRPTS_hrz = length(unique_label);
switch H_hr.Type
case {'sparse'}
Hnum_list_hrz{NRPTS_hrz} = sparse(H_hr.WAN_NUM,H_hr.WAN_NUM);
for i = 1:NRPTS_hrz
vector_list_hrz(i,fin_dir) = unique_dir(i,:);
Hnum_list_hrz{i} = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
end
case {'mat'}
if Num
Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
for i = 1:NRPTS_hrz
vector_list_hrz(i,fin_dir) = unique_dir(i,:);
Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
end
end
if Coe
Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
for i = 1:NRPTS_hrz
vector_list_hrz(i,fin_dir) = unique_dir(i,:);
Hcoe_list_hrz(:,:,i) = sum(Hcoe_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
end
end
end
switch H_hr.Type
case {'sparse'}
H_hr.HnumL = Hnum_list_hrz;
case {'mat'}
if Num
H_hr.HnumL = Hnum_list_hrz;
end
if Coe
H_hr.HcoeL = Hcoe_list_hrz;
end
end
H_hr.vectorL = (vector_list_hrz);
case '2D'
[vector_list_new,sort_label] = sortrows(H_hr.vectorL,fin_dir) ;
[unique_dir,unique_label]= unique(vector_list_new(:,fin_dir),'rows');
switch H_hr.Type
case {'sparse'}
Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label);
case {'mat'}
if Num
Hnum_kf_list_new = Hnum_kf_list(:,:,sort_label);
end
if Coe
Hcoe_kf_list_new = Hcoe_kf_list(:,:,sort_label);
end
end
cutlist(:,1)= unique_label;
cutlist(1:end-1,2)= unique_label(2:end)-1;
cutlist(end,2) = H_hr.NRPTS;
vector_list_hrz = zeros(length(unique_dir),3);
NRPTS_hrz = length(unique_label);
switch H_hr.Type
case {'sparse','num'}
Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
for i = 1:NRPTS_hrz
vector_list_hrz(i,fin_dir) = unique_dir(i,:);
Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
end
case {'mat'}
Hnum_list_hrz = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz);
for i = 1:NRPTS_hrz
vector_list_hrz(i,fin_dir) = unique_dir(i,:);
if Num
Hnum_list_hrz(:,:,i) = sum(Hnum_kf_list_new(:,:,cutlist(i,1):cutlist(i,2)),3);
end
if Coe
Hcoe_list_hrz = sym(zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,NRPTS_hrz));
end
end
end
switch H_hr.Type
case {'sparse'}
H_hr.HnumL = Hnum_list_hrz;
case {'mat'}
if Num
H_hr.HnumL = Hnum_list_hrz;
end
if Coe
H_hr.HcoeL = Hcoe_list_hrz;
end
end
H_hr.vectorL = (vector_list_hrz);
end
end
