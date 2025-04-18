%% useful_tools
clear
useful_matrices(["sigma","tau"]);
%% 
% *变量名不能含有xyz这样的字符*

syms Delta m v lambda u1 u2;
syms k_x k_y k_z c real ;

%% 
% % the coeffs is not allowed use x y z
H_2D = HK(4,2);
H_2D = H_2D ...
    + Term( -Delta + m*k_x^2+m*k_y^2,sigma_z*tau_0)...
    + Term( v*k_x,                   sigma_x*tau_x)...
    + Term( v*k_y,                   sigma_x*tau_y)...
    + Term( lambda*(k_x^2-k_y^2),    sigma_z*tau_x)...
    + Term( 2*lambda*(k_x * k_y),    sigma_z*tau_y)...
    ;

% H_2D = H_2D < 'POSCAR';
H_2D = H_2D.input_Rm();
H_2D = H_2D < 'KPOINTS';
%% 
% % para

b = 1.424e-10;
h_bar =6.582119514e-16; %eV⋅s
v = 7e5 * h_bar /b  ; % m/s
Delta = 0.213 ; % eV
lambda = 0.230; % eV
u1 = 0.26      ;
u2 = 0.23      ;
m = 0.233     ;
c = 3.3  ;
v = 0;
%% 
% 

%v =0;
u1 =0.233;
u2 =0.23;
m =0.19;
Delta =0.213;
tolerance = 0.037;
kmesh = [5,23,5];
H_2Dn = H_2D.Subsall();
%
EIGENCAR = H_2Dn.EIGENCAR_gen();
bandplot(EIGENCAR,[-0.3,0.3]);

%[k_list_cart,klist_f,~]=H_3Dn.findnodes(kmesh,2,tolerance);