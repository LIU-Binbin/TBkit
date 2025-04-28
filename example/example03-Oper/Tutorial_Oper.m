%% Oper Class Tutorial (Oper类教程)
% This tutorial demonstrates basic operations and symmetry element construction using the Oper class
% 本教程演示如何使用Oper类进行基本操作和对称元素构建
% 推荐以mlx打开
%% 1. Constructor Test (构造函数测试)
% 1.1 Create identity operator (创建单位算符)
E = Oper;                     % Default constructor creates identity operator (默认构造函数创建单位算符)
E.pretty('Latex',1);          % Display operator in LaTeX format (以LaTeX格式显示算符)
disp('Identity operator:');   % 显示单位算符
fprintf(E.pretty());          % 文本格式显示

% 1.2 Create inversion operator (创建反演算符)
I = Oper(-eye(3), kron(eye(2), [1,0;0,-1])); % 空间反演 × 自旋操作
disp('Inversion operator:');  % 显示反演算符
fprintf(I.pretty());          % 文本格式显示
I.pretty('Latex',1)         % LaTeX格式显示

% 1.3 Create C3z rotation operator (创建C3z旋转算符)
theta = pi/3;                 % 60度旋转角度
C3z = Oper(Oper.nTheta2Rotation(theta),...      % 生成3D旋转矩阵
          expm(kron([1,0;0,-1],1i*[pi/3,0;0,pi]))); % 自旋旋转操作
disp('C3z rotation operator:'); % 显示C3z算符
fprintf(C3z.pretty());        % 文本格式显示
C3z.pretty('Latex',1)        % LaTeX格式显示

%% 2. Symmetry Operators Construction (对称算符构建)
% 2.1 Time reversal operators (时间反演算符)
tr1 = Oper.time_reversal()           % Default spin-1/2 (默认自旋1/2)
tr2 = Oper.time_reversal(3, nan, 3/2)% 3D空间 + 自旋3/2

% 2.2 Other symmetry operators (其他对称算符)
P = Oper.particle_hole(3,nan)        % 3D粒子空穴对称
Inv = Oper.inversion(3,nan)         % 3D空间反演
Mirror = Oper.mirror([1,0,0])        % 沿x轴镜面反射

% 2.3 Operator sets (算符集合)
A = [tr1, P, Inv];                    % 基础对称操作集合
B = [tr2, P, Inv, Mirror];            % 扩展对称操作集合
% A.generate_group
%% 3. Crystallographic Groups (晶体学群)
% Generate common symmetry groups (生成常见对称群)
Square_group = Oper.square(false,false)     % 二维正方晶系
Cubic_group = Oper.cubic(false,false)       % 三维立方晶系
Hexagonal_2D_group = Oper.hexagonal_2D(false,false) % 二维六角晶系
Hexagonal_group = Oper.hexagonal(false,false)% 三维六角晶系

%% 4. Static Methods (静态方法)
% 4.1 Angle naming convention (角度命名规范)
E = Oper.name_angle(0)                % 0弧度 → "0"
E = Oper.name_angle(0,1)              % 0弧度 (分数π格式)
E = Oper.name_angle(2*pi)             % 自动归一化为0
E = Oper.name_angle(3/2*pi,1)         % 显示为"3/2π"
E = Oper.name_angle(4.7124,1)         % 自动识别为3π/2

% 4.2 Rotation matrix analysis (旋转矩阵解析)
try
    Oper.Rotation2nTheta(1)           % 错误输入测试
catch
    disp('Requires 3×3 rotation matrix'); % 需要3×3矩阵
end
R_matrix = -eye(3);                   % 定义旋转矩阵
[n, theta] = Oper.Rotation2nTheta(R_matrix) % 提取旋转轴和角度

%% 5. Spin and Angular Momentum (自旋与角动量)
% 5.1 Spin matrices generation (自旋矩阵生成)
s = 3/2;                              % 自旋量子数
spin_mat = Oper.spin_matrices(s)      % 生成Sx, Sy, Sz

% 5.2 Orbital angular momentum (轨道角动量矩阵)
L_mat = Oper.L_matrices()             % 默认生成L=1矩阵

% 5.3 Spin rotation operator (自旋旋转算符)
U = Oper.spin_rotation(pi*[0, 1, 0], 3/2) % 绕y轴旋转π (自旋3/2)
disp('Spin rotation matrix:');        % 显示自旋旋转矩阵
disp(U);                              % 输出矩阵数值