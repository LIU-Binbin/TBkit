% 使用该程序对PBANDS进行数据清洗，得到相应轨道在不同自旋方向上的比重
% 该程序目标是为了在自旋极化计算时，打开非线性设置从而规定磁化方向，进一步在能带上投影出自旋上下分量；具体测试并为实现理想功能
% 当前文件夹下需要有：POSCAR，band-EIGENVAL,band-KPOINTS,scf-DOSCAR,band-PBAND_XX_SOC_XX.dat

filename1 = "PBAND_Mn_SOC_SY.dat";
[WEIGHTCAR_3d_mat1,~,~] = WEIGHTCAR_gen(filename1,-1,'vaspkit-band-silence');
Name_list1 = "PBAND_Mn_SOC_SY"; % 这个命名无所谓
width1= size(WEIGHTCAR_3d_mat1,3);
WEIGHTCAR_struct1= temp_3d_mat2struct(WEIGHTCAR_3d_mat1,Name_list1,width1);

% filename2 = "PBAND_Mn_SOC_SY.dat";
% [WEIGHTCAR_3d_mat2,~,~] = WEIGHTCAR_gen(filename2,-1,'vaspkit-band-silence');
% Name_list2 = "PBAND_Mn_SOC_SY"; % 这个命名无所谓
% width2= size(WEIGHTCAR_3d_mat2,3);
% WEIGHTCAR_struct2= temp_3d_mat2struct(WEIGHTCAR_3d_mat2,Name_list2,width2);


% % WEIGHTCAR = WEIGHTCAR_struct1(5).WEIGHTCAR() + WEIGHTCAR_struct1(6).WEIGHTCAR() + WEIGHTCAR_struct1(7).WEIGHTCAR()...
%    + WEIGHTCAR_struct1(8).WEIGHTCAR() + WEIGHTCAR_struct1(9).WEIGHTCAR()...
%    +WEIGHTCAR_struct2(5).WEIGHTCAR() + WEIGHTCAR_struct2(6).WEIGHTCAR() + WEIGHTCAR_struct2(7).WEIGHTCAR()...
%    + WEIGHTCAR_struct2(8).WEIGHTCAR() + WEIGHTCAR_struct2(9).WEIGHTCAR();
% WEIGHTCAR = WEIGHTCAR_struct1(10).WEIGHTCAR() +WEIGHTCAR_struct2(10).WEIGHTCAR();

%WEIGHTCAR = WEIGHTCAR_struct1(10).WEIGHTCAR();

 WEIGHTCAR = WEIGHTCAR_struct1(5).WEIGHTCAR() + WEIGHTCAR_struct1(6).WEIGHTCAR() + WEIGHTCAR_struct1(7).WEIGHTCAR()...
    + WEIGHTCAR_struct1(8).WEIGHTCAR() + WEIGHTCAR_struct1(9).WEIGHTCAR();
EIGENCAR = EIGENVAL_read(); % 没有改为新计算的生成文件
pbandplot(WEIGHTCAR,EIGENCAR) %这里两组数据如果维度不一致，似乎不会报错？
ylim([-3,1]);

function WEIGHTCAR_struct = temp_3d_mat2struct(WEIGHTCAR_3d,Name,width)
import vasplib_tool.*
% f or not
if width >10
    num2orbital_name =  orbital_maprule_gen(1);
else
    num2orbital_name = orbital_maprule_gen(0);
end
for i = 1:width
    WEIGHTCAR_struct(i).WEIGHTCAR = WEIGHTCAR_3d(:,:,i);
    WEIGHTCAR_struct(i).displayname = Name + "-"+num2orbital_name(i);
end
end