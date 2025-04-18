%% 拓扑半金属图示
% an illustration of the drumhead state in nodal-line semi-metals as compared 
% to related topological semi-metal counterparts.
% 
% 
%% WSM
% wsm 直接采用 QWZ模型
% 
% 

s_0   = pauli_matrix(0)  ;  s_x = pauli_matrix(1);  s_y =  pauli_matrix(2) ;  s_z = pauli_matrix(3);
sigma_0 = pauli_matrix(0);sigma_x =  pauli_matrix(1);sigma_y =  pauli_matrix(2);sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0);  tau_x =  pauli_matrix(1);  tau_y =  pauli_matrix(2);  tau_z = pauli_matrix(3);

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A0 real;
syms k_x k_y k_z real;
% kp

M       = M0-M2*(k_x^2+k_y^2)-M1*(k_z^2);
E0k     = C0+C2*(k_x^2+k_y^2)+C1*(k_z^2);
A       = A0;%+A2*(k_x^2+k_y^2)+A1*(k_z^2)
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
QWZ_3D = HK(2,2);
QWZ_3D = QWZ_3D ...
    +Term(A*k_x ,tau_x )...
    +Term(A*k_y ,tau_y )...
    +Term(E0k   ,tau_0 )...
    +Term(M     ,tau_z );...
QWZ_3D = QWZ_3D.Subsall('sym');
QWZ_3D = QWZ_3D <'POSCAR_2';
%%
QWZ_3D_kp = QWZ_3D.sym()
%QWZ_3D_kp_pauli = QWZ_3D.pauliDecomposition
% kp2TB

QWZ_3D_TB= QWZ_3D.kp2TB();
QWZ_3D_TB.list();
% Para

% unit eV 
M0     = 1   ;
%M0     = -0.25   ;
% unit eV Ang
A0      =  1      ;
%B      =  4.1      ;
% unit eV Ang^2
C0     =  0;
C1     =  0.1     ;
C2     =  0     ;
M1     =  0.5      ;
M2     =  0.5  ;
% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% NumTB

QWZ_3D_TB_n = QWZ_3D_TB.Subsall();
QWZ_3D_TB_n = QWZ_3D_TB_n <'KPOINTS_4';
%% Bulk
% bandstructrue

EIGENCAR = QWZ_3D_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-2,2],'title',"WSM-TETRA",'POSCAR','POSCAR_4','KPOINTS','KPOINTS_4');
% 3d Band

[EIGENCAR_3D,klist1,klist2] = QWZ_3D_TB_n.EIGENCAR_gen_3D([101,101],[-0.5,-0.0,-0.5;1,0,0;0 0 1],'fin_dir',2);
% [fig,ax] = QWZ_3D_TB_n.BZplot(QWZ_3D_TB_n.Rm,'mode','3D');
%%
[~,Ax]= Figs(1,3);
bandplot3d(EIGENCAR_3D(:,:,1),klist1,klist2,'ax',Ax(1),'WEIGHTCAR',-ones(101),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
bandplot3d(EIGENCAR_3D(:,:,2),klist1,klist2,'ax',Ax(1),'WEIGHTCAR',ones(101),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
%% DSM
% Dsm 单点直接采用 
% 
% 

syms A k_x k_y real;

%
DSM_3D = HK(4,2);
DSM_3D = DSM_3D ...
    +Term(A*k_x ,sigma_z*tau_x )...
    +Term(A*k_y ,sigma_0*tau_y );
DSM_3D = DSM_3D.Subsall('sym');
DSM_3D = DSM_3D <'POSCAR_4';
%%
DSM_3D_kp = DSM_3D.sym()

DSM_3D_TB= DSM_3D.kp2TB();

% Para

% unit eV 
A     = 2  ;

% unit    Ang
a      =  1        ;
b      =  1       ;
c      =  1        ;
% NumTB

DSM_3D_TB_n = DSM_3D_TB.Subsall();
DSM_3D_TB_n = DSM_3D_TB_n <'KPOINTS_4';
%% Bulk
% bandstructrue

EIGENCAR = DSM_3D_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-2,2],'title',"DSM-TETRA",'POSCAR','POSCAR_4','KPOINTS','KPOINTS_4');
% 3d Band

[EIGENCAR_3D_DSM,klist1,klist2] = DSM_3D_TB_n.EIGENCAR_gen_3D([101,101],[-0.25,-0.25,-0.0;0.5,0,0;0 0.5 0],'fin_dir',3);
%%


bandplot3d(EIGENCAR_3D_DSM(:,:,1)-0.1,klist1,klist2,'ax',Ax(2),'WEIGHTCAR',-ones(101)-1,'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
bandplot3d(EIGENCAR_3D_DSM(:,:,2),klist1,klist2,'ax',Ax(2),'WEIGHTCAR',-ones(101),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
bandplot3d(EIGENCAR_3D_DSM(:,:,3),klist1,klist2,'ax',Ax(2),'WEIGHTCAR',ones(101),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
bandplot3d(EIGENCAR_3D_DSM(:,:,4)+0.1,klist1,klist2,'ax',Ax(2),'WEIGHTCAR',ones(101)+1,'xlabel','k_x','ylabel','k_y', ...
    'cmap',[0,0,1;0.2,0,0.6;0.6,0,0.2;1,0,0],'FaceAlpha',0.5);

%% 
% NLSM
% 
% 
% 
% 10.1088/1674-1056/25/11/117106

s_0   = pauli_matrix(0)  ;  s_x = pauli_matrix(1);  s_y =  pauli_matrix(2) ;  s_z = pauli_matrix(3);
sigma_0 = pauli_matrix(0);sigma_x =  pauli_matrix(1);sigma_y =  pauli_matrix(2);sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0);  tau_x =  pauli_matrix(1);  tau_y =  pauli_matrix(2);  tau_z = pauli_matrix(3);

syms C0 C1 C2 real;
syms M0 M1 M2 real;
syms A0 A real;
syms k_x k_y k_z real;
% kp

M       = M0*(1-(k_x^2+k_y^2));
k_plus  = k_x + 1i* k_y;
k_minus = k_x - 1i* k_y;
%
NLSM_3D = HK(2,2);
NLSM_3D = NLSM_3D ...
    +Term(A*k_z ,tau_x )...
    +Term(M     ,tau_z );...
NLSM_3D = NLSM_3D.Subsall('sym');
NLSM_3D = NLSM_3D <'POSCAR_2';
%%
NLSM_3D_kp = NLSM_3D.sym()
%NLSM_3D_kp_pauli = NLSM_3D.pauliDecomposition
% kp2TB

NLSM_3D_TB= NLSM_3D.kp2TB();


Hsym = rewrite(NLSM_3D_TB.sym("simple",1),'sincos');
simplify(eig(subs(Hsym,k_z,0)))
% NLSM_3D_TB.list();
% Para

% unit eV 
M0     = 0.5   ;
%M0     = -0.25   ;
% unit eV Ang
A      =  1     ;
% unit    Ang
a      =  3        ;
b      =  3       ;
c      =  3        ;
% NumTB

NLSM_3D_TB_n = NLSM_3D_TB.Subsall();
NLSM_3D_TB_n = NLSM_3D_TB_n <'KPOINTS_4';
%% Bulk
% bandstructrue

EIGENCAR = NLSM_3D_TB_n.EIGENCAR_gen();
bandplot(EIGENCAR,[-2,2],'title',"nlSM-TETRA",'POSCAR','POSCAR_4','KPOINTS','KPOINTS_4');
% 3d Band

[EIGENCAR_3D,klist1,klist2] = NLSM_3D_TB_n.EIGENCAR_gen_3D([301,301],[-0.5,-0.5,-0.0;1,0,0;0 1 0],'fin_dir',3);
% [fig,ax] = NLSM_3D_TB_n.BZplot(NLSM_3D_TB_n.Rm,'mode','3D');

%%

bandplot3d(EIGENCAR_3D(:,:,1),klist1,klist2,'ax',Ax(3),'WEIGHTCAR',-ones(301),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
bandplot3d(EIGENCAR_3D(:,:,2),klist1,klist2,'ax',Ax(3),'WEIGHTCAR',ones(301),'xlabel','k_x','ylabel','k_y','cmap',[0,0,1;1,0,0],'FaceAlpha',0.5);
view(Ax(1),105,12);
axis(Ax(1),'off');
view(Ax(2),105,12);
axis(Ax(2),'off');
view(Ax(3),105,12);
axis(Ax(3),'off');
%%
export_fig(gcf,'temp.png','-png','-r600');
%% 
% 
% 
%