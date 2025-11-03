% 重复 test
% 2025.10.28，Yi-Qun You
clear;
clc;
%% 
s_0   = pauli_matrice(0)  ;  s_x = pauli_matrice(1);  s_y =  pauli_matrice(2) ;  s_z = pauli_matrice(3);
sigma_0 = pauli_matrice(0);sigma_x =  pauli_matrice(1);sigma_y =  pauli_matrice(2);sigma_z = pauli_matrice(3);
tau_0   = pauli_matrice(0);  tau_x =  pauli_matrice(1);  tau_y =  pauli_matrice(2);  tau_z = pauli_matrice(3);

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A0 B real;
syms k_x k_y k_z real;

%% 
M       = -M0+M2*(k_x^2+k_y^2)+M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
A       = A0;
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;

BHZ = HK(4,2); % 4带，展开到k^2
BHZ = BHZ ...
    +Term(A*k_x ,kron(sigma_z,tau_x) )...
    +Term(A*k_y ,kron(sigma_0,tau_y) )...
    +Term(E0k   ,kron(sigma_0,tau_0 ))...
    +Term(M     ,kron(sigma_0,tau_z) )...
    +Term(B*(k_x^2-k_y^2),kron(sigma_x,tau_x) )...
    +Term(2*B*(k_x*k_y),kron(sigma_y,tau_x ))...
    ;
sym(BHZ)
%% 
M0 = -1;
A0 =  1;
C0 =  0; 
C1 =  0;
C2 =  0;
M1 = -0.5;
M2 = -0.5;
B  =  1;

% a = 1; 
% b = 1;
% c = 1;

BHZ = BHZ <'POSCAR_6'; % KP的POSCAR信息并不重要
%% 
BHZ_TB= BHZ.kp2TB(); % 没有加其他转化的要求
BHZ_TB.list('vpa',false); % 考虑HK2Htrig2HR，具体矩阵形式约束最终为9个跃迁方向
%% 
BHZ_Htrig = BHZ_TB.HR2Htrig();
sym(BHZ_Htrig)
%% 
BHZ_TB_n = BHZ_TB.Subsall();
BHZ_TB_n = BHZ_TB_n <'KPOINTS_6';
EIGENCAR = BHZ_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-2,2],'title',"HODSM-HEX",'KPOINTS','KPOINTS_6'); % HEX表示六角晶格
%% 
[BFCAR,~,klist_l] = BHZ_TB_n.WilsonLoop();
[fig,ax] = WilsonLoopPlot(BFCAR,klist_l); 