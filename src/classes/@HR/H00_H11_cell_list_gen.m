function [H00_H11_cell_list_1,H00_H11_cell_list_2] = H00_H11_cell_list_gen(H_hr,fin_dir,principle_layer)
if nargin <3
principle_layer = 3;
end
if nargin <2
fin_dir = 2;
end
k_n=length(H_hr.klist_frac(:,1));
H00_H11_cell_list_1 = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,k_n);
H00_H11_cell_list_2 = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,k_n);
pb = TBkit_tool_outer.CmdLineProgressBar('Hamiltionian Generating ');
for j = 1:k_n
pb.print(j,k_n);
kpoints_f = H_hr.klist_frac(j,:);
H_hrz = H_hr.hrz_gen(kpoints_f,fin_dir);
[H00,H01,~] = H_hrz.Green_prepare(principle_layer,fin_dir);
H00_H11_cell_list_1(:,:,j) = H00;
H00_H11_cell_list_2(:,:,j) = H01;
end
pb.delete();
end
