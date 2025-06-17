%% BHZ model
% Useful matrix
%% useful tool

s_0   = pauli_matrix(0)  ;  s_x = pauli_matrix(1);  s_y =  pauli_matrix(2) ;  s_z = pauli_matrix(3);
sigma_0 = pauli_matrix(0);sigma_x =  pauli_matrix(1);sigma_y =  pauli_matrix(2);sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0);  tau_x =  pauli_matrix(1);  tau_y =  pauli_matrix(2);  tau_z = pauli_matrix(3);
% k·p model

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A B real;
syms k_x k_y k_z real;
%
M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
BHZ = HK(4,3)
%
BHZ = BHZ ...
    +Term(A*k_x ,sigma_z*tau_x )...
    +Term(A*k_y ,sigma_0*tau_y )...
    +Term(E0k   ,sigma_0*tau_0 )...
    +Term(M     ,sigma_0*tau_z )...
    +Term(B*k_z*(k_x^2-k_y^2),sigma_x*tau_x )...
    +Term(2*k_z*B*(k_x*k_y),sigma_y*tau_x )...
    ;
%
BHZ_PT = HK(4,2)
%
BHZ_PT = BHZ_PT ...
    +Term(A*k_x ,sigma_z*tau_x )...
    +Term(A*k_y ,sigma_0*tau_y )...
    +Term(E0k   ,sigma_0*tau_0 )...
    +Term(M     ,sigma_0*tau_z )...
    +Term(B*(k_x^2-k_y^2),sigma_x*tau_x )...
    +Term(2*B*(k_x*k_y),sigma_y*tau_x )...
    ;
% Oper

Cubic = Oper.cubic(false,false);
C4 = Cubic(17);
C4.U = expm(diag([1i*pi/4,1i*3*pi/4,-1i*pi/4,-1i*3*pi/4]));
I = Oper.inversion(3,diag([1,-1,1,-1]));
%%
BHZ.applyOper(C4)
%%
BHZ.applyOper(I)
%%
BHZ_PT
BHZ_PT.symvar_list
%%
BHZ_PT_after_P = BHZ_PT.applyOper(I)
BHZ_PT_after_P.symvar_list

%% 
% 

BHZ_building1 = HK(4,2);
BHZ_building1 = BHZ_building1 <= 'POSCAR';
J_mat = BHZ_building1.Jmat_gen();

Tr = Oper.time_reversal(3,double(-1i*sigma_y*tau_0),nan)
My = Oper.mirror([0,1,0],double(1i*sigma_y*tau_z),nan)
I = Oper.inversion(3,diag([1,-1,1,-1]))

C4z = Oper.rotation(1/4,[0,0,1],false,sym(expm(1i*pi/4 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true)
C6z = Oper.rotation(1/6,[0,0,1],false,sym(expm(1i*sym(pi)/6 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true)
C3z = Oper.rotation(1/3,[0,0,1],false,sym(expm(1i*sym(pi)/3 *double(sigma_z*(tau_0-2*tau_z)))),nan,'sym',true)
groups = [Tr,My,I,C4z];
groups2 = [Tr,My,I,C6z];
groups3 = [Tr,My,I,C3z];
BHZ_building_C4 = BHZ_building1.applyOper(groups)
BHZ_building_C6 = BHZ_building1.applyOper(groups2)
BHZ_building_C6 = BHZ_building1.applyOper(groups2)
BHZ_building_C6 = BHZ_building1.applyOper(groups2)
BHZ_building_C3 = BHZ_building1.applyOper(groups3)

syms k_x k_y k_z k_plus k_minus real;
%% 
%