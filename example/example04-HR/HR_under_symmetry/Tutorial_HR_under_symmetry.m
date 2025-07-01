%% 对称性约束下的紧束缚模型构建 | Tight-binding Model Construction under Symmetry Constraints
%% 理论背景 | Theoretical Background
% 对称操作下哈密顿量的变换公式 | Hamiltonian transformation under symmetry operations $$\widetilde{H}^{i 
% j}\left[\mathbf{R}^{\prime}\right]=\frac{1}{|G|} \sum_{{g  \in G }}^{lm} D_{i 
% l}(g) H^{l m}\left[S_g^{-1}\left(\mathbf{R}^{\prime}-\mathbf{T}_{i  j}^{m l}\right)\right] 
% D_{m j}\left(g^{-1}\right)$$
% 
% where $S_g$ is the real-space rotation matrix of the symmetry $g, \mathbf{T}_{i    
% j}^{m l}=S_g\left(\mathbf{t}_m-\mathbf{t}_l\right)-\mathbf{t}_j-\mathbf{t}_i$, 
% and the indices $m, l$ only go over values for which $\mathbf{T}_{i j}^{m l}$ 
% is a lattice vector. 

% 参数说明 | Parameters:
% - S_g: 实空间旋转矩阵 | Real-space rotation matrix
% - T_{ij}^{ml} = S_g(t_m - t_l) - t_j - t_i: 位置修正项 | Position correction term
% - m,l: 仅当T_{ij}^{ml}为晶格矢量时有效 | Valid only when T_{ij}^{ml} is lattice vector
%% 对称性生成紧束缚模型 | Symmetry-Generated TB Model
% 初始化HR对象 | Initialize HR object

tic;
Graphene = HR(2);  % 二维系统 | 2D system
Graphene = 'POSCAR' > Graphene;  % 读取晶体结构 | Read crystal structure

% 最近邻搜索 | Nearest neighbor search
Graphene = Graphene.nn([1,1,0], 4, 1.15);  % 搜索参数：方向容差/截断半径 | Search parameters: direction tolerance/cutoff

% 模型初始化 | Model initialization
Graphene = Graphene.init('level_cut', 1, 'onsite', 1);
disp('初始哈密顿量 | Initial H:');
list(Graphene);
Graphene.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);

% 定义对称操作 | Define symmetry operations
C3 = Oper.rotation(1/3, [0,0,1], false, eye(2));  % C3旋转
Tr = Oper.time_reversal(3, diag([1, 1]));        % 时间反演对称性
Chiral = Oper.chiral(3, diag([1, -1]));          % 手性对称性
Groups = generate_group([C3, Tr, Chiral]);       % 生成对称群

% 应用对称约束 | Apply symmetry constraints
Graphene_test = Graphene.applyOper(Groups, 'generator', true);
disp('对称约束后的哈密顿量 | H after Symmetry Constraints:');
list(Graphene_test);  % 显示非零hopping项 | Display non-zero hopping terms

% 可视化结果 | Visualization
Graphene_test.show('HOPPING', 'scale', 2.456, 'atomscale', 0.1, 'TwoD', true);

time = toc;
fprintf('Time cost %f s\n',time);
% other info
[Graphene_test2,Sublist,Unique_term] = Graphene_test.unique();
Graphene_test2.sym
%
Graphene_test3 = Graphene.applyOper(C3);
Graphene_test3.list;
Graphene_test3.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
Graphene_test3 = Graphene_test3.applyOper(C3,'generator',1);
Graphene_test3.list;
Graphene_test3.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
Graphene_test3 = Graphene_test3.applyOper(Tr,'generator',1);
Graphene_test3.list;
%% 快速参数化模式 | Fast Parameterization Mode
% 初始化快速模式 | Initialize fast mode

clear;
tic;
Graphene_fast = HR(2);
Graphene_fast =  'POSCAR' > Graphene_fast;

% 最近邻搜索 | Efficient nearest neighbor search
Graphene_fast = Graphene_fast.nn([1,1,0], 1e-4, 1.15);

% 快速参数初始化 | Fast parameter initialization
Graphene_fast = Graphene_fast.init('level_cut', 1, 'onsite', 1, 'fast', true);

% restore in CvectorL
Graphene_fast.list();
%
C3 = Oper.rotation(1/3,[0,0,1],false,eye(2))
Tr = Oper.time_reversal(3,diag([1 1]))
Chiral = Oper.chiral(3,diag([1,-1]))
%
Groups = generate_group([C3,Tr,Chiral])
% 快速应用对称操作 | Fast symmetry application
Graphene_fast_sym = Graphene_fast.applyOper(Groups, 'generator', true, 'fast', true);

% 生成hop | Generate 
Graphene_fast_vis = Graphene_fast_sym.GenfromOrth(); 

% 可视化结果 | Visualization
Graphene_fast_vis.show('HOPPING', 'scale', 2.456, 'atomscale', 0.1, 'TwoD', true);

time = toc;
fprintf('Faster Time cost %f s\n',time);
%% Kane-Mele模型案例 | Kane-Mele Model Case
% 初始化Kane-Mele模型 | Initialize Kane-Mele model

clear;
KaneMele = HR(4);
KaneMele = KaneMele < 'POSCAR_KM';  % 读取拓扑绝缘体结构 | Read topological insulator structure
KaneMele.Rm = sym(KaneMele.Rm);
KaneMele.orbL = sym(KaneMele.orbL);
KaneMele= KaneMele.nn([1,1,0],1e-2,1.15);
%
KaneMele = KaneMele.init('level_cut',2,"onsite",1,'fast',true);

% 定义复杂对称操作 | Define complex symmetries
C3 = Oper.rotation(1/3, [0,0,1], false, double(expm(-1i*(pi/3)*gamma_matrix(2,4)))); % 自旋轨道耦合旋转
Tr = Oper.time_reversal(3, double(-1i*gamma_matrix(4,5)));  % 时间反演
I = Oper.inversion(3, double(-gamma_matrix(1)));            % 空间反演
Mx = Oper.mirror([1,0,0], double(1i*gamma_matrix(2,5)));     % x镜像
My = Oper.mirror([0,1,0], double(1i*gamma_matrix(2,3)));     % y镜像

% 应用对称约束 | Apply symmetry constraints
Groups_KM = generate_group([C3, Tr, I, Mx, My]);
KaneMele_test = KaneMele.applyOper([C3,I,Mx,My,Tr], 'generator', true, 'fast', true);
KaneMele_test = KaneMele_test.GenfromOrth('Accuracy',1e-6);
% 参数替换与化简 | Parameter substitution and simplification
syms t lambda_SO E_pz T_2 real;
T_2 = 0; % omit real component <<i,j>>
Varlist = KaneMele_test.symvar_list;
KaneMele_test2 = subs(KaneMele_test,Varlist,[T_2,t,E_pz,lambda_SO]);
disp('最终化简形式 | Final Simplified Form:');
disp(simplify(rewrite(KaneMele_test2.sym(), 'sincos')));
%% 
t  = 1;
lambda_SO = 0.06;
E_pz = 0;
T_2 = 1;
[~,Ax] = Figs(2,2);
KaneMele_test2_n = KaneMele_test2.Subsall();
EIGENCAR_origin = KaneMele_test2_n.EIGENCAR_gen();

bandplot(EIGENCAR_origin,'ax',Ax(1),'title','KM:t = 1, \lambda_{SO} = 0.06');
KaneMele_test2_n.HnumL = KaneMele_test2_n.HnumL + 0.1*rand(size(KaneMele_test2_n.HnumL));


KaneMele_test2_n.bandplot('ax',Ax(2),'title','KM:t = 1, \lambda_{SO} = 0.06,Error:0.1')
KaneMele_test3_n = KaneMele_test2_n.applyOper([C3,I,Mx,My,Tr], 'generator', true);
EIGENCAR_sym = KaneMele_test3_n.EIGENCAR_gen();
bandplot(EIGENCAR_sym,'ax',Ax(3),'title','KM:t = 1, \lambda_{SO} = 0.06,Error:0.1,sym');
bandplot({EIGENCAR_origin,EIGENCAR_sym},'ax',Ax(4),'title','KM-vs-sym','Color',[1 0 0;0 0 1],'legends',["KM","KM-sym"]);
for i =1:4
    axis(Ax(i),'normal');
end
