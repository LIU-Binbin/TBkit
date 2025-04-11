function H_hr = tbbox_in_gen(H_hr,options)
% TBBOX_IN_GEN Generate input file for tight-binding calculations
%
%   H_HR = TBBOX_IN_GEN(H_HR,OPTIONS) creates tbbox.in input file
%
%   Inputs:
%       H_hr - HR object
%       options.Accuracy - Rounding accuracy [default: 1e-6]
%       options.OperL - Operator list [default: []]
%       options.hr - Generate hr file [default: true]
%       options.wan - Wannier mode [default: false]
%       options.lda - LDA mode [default: true]
%   Output:
%       H_hr - Processed HR object
%
%   Notes:
%       - Creates tbbox.in and hr.dat files
%       - Handles SOC and non-SOC cases
%       - Uses spglib for symmetry detection
arguments
H_hr HR;
options.Accuracy = 1e-6;
options.OperL = [];
options.hr = true;
options.wan = false;
options.lda = true;
end
import spglib_matlab.*;
Accuracy = -log10(options.Accuracy);
%fprintf('Reseq the orbital by your self');
[unique_orbL,unique_seqL,unique_seqL_inv] = unique(H_hr.orbL,'stable','rows');
unique_orbL_counts = accumarray(unique_seqL_inv,1);
lattice = H_hr.Rm;
position = unique_orbL.';
types = H_hr.quantumL(unique_seqL,1);
if isempty(options.OperL )
try
SpglibDataset  = spglib_matlab.spg_get_dataset_from_sites(lattice,position,types);
rotations = double(SpglibDataset.rotations);
translations = double(SpglibDataset.translations);
n_operations = SpglibDataset.n_operations;
catch
rotations = eye(3);
translations = zeros(1,3);
n_operations = 1;
end
else
n_operations = length(options.OperL);
rotations = zeros(3,3,n_operations);
translations  = zeros(1,3,n_operations);
for i = 1:n_operations
rotations(:,:,i) = options.OperL(i).R;
translations(:,:,i) = options.OperL(i).t;
end
end
sym_oper_list = zeros(n_operations,9);
for i = 1:n_operations
sym_oper_list(i,1) = i;
sym_oper_list(i,2) = det(rotations(:,:,i));
[n,theta]= spglib_matlab.Rotation2nTheta(rotations(:,:,i),H_hr.Rm);
sym_oper_list(i,3) = roundn(theta,-Accuracy);
sym_oper_list(i,4) = n(1);
sym_oper_list(i,5) = n(2);
sym_oper_list(i,6) = n(3);
sym_oper_list(i,7) =roundn(translations(1,1,i),-Accuracy);
sym_oper_list(i,8) =roundn(translations(1,2,i),-Accuracy);
sym_oper_list(i,9) =roundn(translations(1,3,i),-Accuracy);
end
if sum(H_hr.quantumL(:,4)) > 0 || options.lda
casename  = 'lda';
SOCcounts = 1;
else
casename  = 'soc';
SOCcounts = 2;
end
filename_hr = casename+"_hr.dat";
ntau = length(unique_seqL);
fileID = fopen('tbbox.in','w');
fprintf(fileID,' case = %s !lda or soc ! spinless of spinful ！dont use soc or consider\n\n',casename);
fprintf(fileID,' proj:\n');
%fprintf(fileID,'orbt = 1  !Orbital convertion 1: s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2 !!!\n');
fprintf(fileID,' orbt = 1\n');
%fprintf(fileID,'ntau = %d !Orbital convertion 2: s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy (Wannier90) !!!\n',ntau);
fprintf(fileID,' ntau = %d\n',ntau);
for i  = 1:ntau
fprintf(fileID,' %11.8f %11.8f %11.8f %2d %2d',...
unique_orbL(i,1),unique_orbL(i,2),unique_orbL(i,3),...
H_hr.quantumL(unique_seqL(i),1),unique_orbL_counts(i)/SOCcounts);
switch i
case 1
fprintf(fileID,' !! x1, x2, x3, itau, iorbit \n');
case 2
fprintf(fileID,' !! ! itau 指标第几个原子,iorbit 指标什么类型的轨道 \n');
case 3
fprintf(fileID,' ! 1 s (ir2tb_pz 1 pz)  \n');
case 4
fprintf(fileID,' ! 2 (ir2tb_pz 2 s pz)  \n');
case 5
fprintf(fileID,' ! 3 px py pz  \n');
case 6
fprintf(fileID,' ! 4 s px py pz  \n');
case 7
fprintf(fileID,' ! 5 d  \n');
case 8
fprintf(fileID,' ! 6 s d \n');
case 9
fprintf(fileID,' ! 7 f  \n');
case 10
fprintf(fileID,' ! 8 p d  \n');
case 11
fprintf(fileID,' ! 9 s p d  \n');
otherwise
fprintf(fileID,'\n');
end
end
fprintf(fileID,' end projections\n\n');
fprintf(fileID,' kpoint:\n');
fprintf(fileID,' kmesh = 1\n');
fprintf(fileID,' Nk = 8  !8 lines\n');
for i =1:8
switch i
case 1
b1 =0;b2 =0;b3 =0;
case 2
b1 =0.5;b2 =0;b3 =0;
case 3
b1 =0;b2 =0.5;b3 =0;
case 4
b1 =0.5;b2 =0.5;b3 =0;
case 5
b1 =0;b2 =0;b3 =0.5;
case 6
b1 =0.5;b2 =0;b3 =0.5;
case 7
b1 =0;b2 =0.5;b3 =0.5;
case 8
b1 =0.5;b2 =0.5;b3 =0.5;
end
fprintf(fileID,' %11.8f %11.8f %11.8f !k%d :b1 b2 b3 \n',b1,b2,b3,i);
end
fprintf(fileID,' end kpoint_path\n\n');
Gk = eye(3)/H_hr.Rm;
fprintf(fileID,' unit_cell:\n');
for i =H_hr.Dim
fprintf(fileID,'    %12.9f %12.9f %12.9f    %12.9f %12.9f %12.9f \n',...
H_hr.Rm(i,1),H_hr.Rm(i,2),H_hr.Rm(i,3),Gk(1,i),Gk(2,i),Gk(3,i));
end
for i =1 :size(sym_oper_list,1)
fprintf(fileID,' %4d ',sym_oper_list(i,1));
%fprintf('%2d ',sym_oper_list(i,1));
for j =2:9
fprintf(fileID,' %11.8f ',sym_oper_list(i,j));
%fprintf('%11.7f ',sym_oper_list(i,j));
end
fprintf(fileID,'\n');
end
fprintf(fileID,' end unit_cell_cart\n');
if options.hr
if options.wan
Norb = H_hr.WAN_NUM;
H_hr = H_hr.reseq([(1:Norb/2)*2-1,(1:Norb/2)*2]);
end
H_hr.Gen_hr(filename_hr);
end
fclose(fileID);
end
