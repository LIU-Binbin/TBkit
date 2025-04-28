%% Group Theory Tutorial/群论教程
%% Group Class Test/群类测试
%% Basic Concepts of Groups/群的基本概念
% Group Definition/群的定义:
%% 
% ✦ A group is a non-empty set with a binary operation (called group multiplication) that satisfies the following conditions:
% 群是定义了二元运算（称为群乘法）且满足下列条件的非空集合：
%
% (1). Closure: For all a,b ∈ G, ab ∈ G.
% 封闭性：对任意a,b ∈ G，满足ab ∈ G
% 
% (2). Associativity: For all a,b,c ∈ G, (ab)c = a(bc)
% 结合律：对任意a,b,c ∈ G，满足(ab)c = a(bc)
% 
% (3). Identity Element: There exists a unique e ∈ G such that ea = ae = a for all a ∈ G
% 存在唯一单位元e ∈ G，使得对所有a ∈ G有ea = ae = a
% 
% (4). Inverse Element: For each a ∈ G, there exists a unique b ∈ G such that ab = ba = e
% 对每个a ∈ G，存在唯一逆元b ∈ G使得ab = ba = e

% Create Pauli matrices group/生成泡利矩阵群
sigma_x = group(([0 1;1 0]));  % Pauli X matrix/泡利X矩阵
sigma_y = group(([0 1i;-1i 0])); % Pauli Y matrix/泡利Y矩阵
sigma_z = group(([1 0;0 -1]));  % Pauli Z matrix/泡利Z矩阵
Pauli_class = [sigma_x,sigma_y,sigma_z]; % Initial group elements/初始群元素

%% Closure Check/封闭性检查
disp('Original group closure status/原始群封闭状态:')
ismember(sigma_x,Pauli_class)
ismember(sigma_y,Pauli_class)
ismember(sigma_x*sigma_y,Pauli_class) % Should be false/应为假
Pauli_class.closure % Verify closure property/验证封闭性

% Generate closed group/生成封闭群
Pauli_class = Pauli_class.generate_group(); % Generate complete group/生成完整群
disp('After group generation/生成群之后:')
Pauli_class.closure % Now should be closed/现在应为封闭

%% Associativity Check/结合律验证
check_associative = sigma_x*(sigma_y*sigma_z) == (sigma_x*sigma_y)*sigma_z;
disp(['Associativity holds: ' num2str(check_associative)])

%% Identity Element Check/单位元验证
identity_check = Pauli_class(1).E == Pauli_class; % E is identity element/E是单位元
disp('All elements equal to identity when multiplied by E:')
disp(identity_check)

%% Inverse Element Check/逆元验证
inv_check = Pauli_class(2).inv == Pauli_class; % Find inverse elements/寻找逆元
disp('Inverse element check:')
disp(inv_check)

inv_operation_check = Pauli_class(2).inv * Pauli_class(2) == Pauli_class(1).E;
disp(['Inverse operation valid: ' num2str(inv_operation_check)])

%% Group Order/群的阶数
disp(['Group order: ' num2str(Pauli_class.order)])

%% Rearrangement Theorem/重排定理
% 定理内容：对于群G中的任意元素g，左乘（或右乘）g后的集合g*G = {g*h | h∈G} 仍是原群的排列。
% Theorem: For any element g in group G, the set g*G = {g*h | h∈G} is a permutation of G.

rearrangement_check = Pauli_class(randi([1,16]))*Pauli_class == Pauli_class;
disp(['Rearrangement theorem holds: ' num2str(all(rearrangement_check))])

%% Subgroups/子群
% Every group has at least two subgroups/每个群至少有两个子群
disp('Trivial subgroup checks:')
all(ismember(Pauli_class(1).E, Pauli_class)) % Trivial subgroup/平凡子群
all(ismember(Pauli_class, Pauli_class))      % Improper subgroup/非真子群

% Generate all subgroups/生成所有子群
SubClassOfPauli_class = Pauli_class.subgroup();
disp(['Number of subgroups found: ' num2str(length(SubClassOfPauli_class))])

%% Cosets/陪集
H = SubClassOfPauli_class{6};       % Select a subgroup/选择子群
gH = Pauli_class(2)*H;              % Left coset/左陪集
Hg = H*Pauli_class(2);              % Right coset/右陪集

% Coset theorem verification/陪集定理验证
gH2 = Pauli_class(3)*H;
disp(['Cosets are either equal or disjoint: ' num2str(gH2 == gH)])

%% Conjugation Classes/共轭类
Square_group = Oper.square(false,false); % Create C4v group/创建C4v群
G_alpha = Square_group(5);           % Select group element/选择群元素
G_beta = Square_group(6);

% Check conjugation/检查共轭
conj_check = conjugate(G_alpha, G_beta, Square_group);
disp(['Elements are conjugate: ' num2str(conj_check)])

% Get conjugation class/获取共轭类
G_alphaConjugation = G_alpha.conjugation(Square_group);
disp('Conjugation class elements:')
disp(G_alphaConjugation)

%% Normal Subgroups/正规子群
NormalSubgroup = normalsubgroup(Square_group); % Find normal subgroups/寻找正规子群
disp(['Number of normal subgroups: ' num2str(length(NormalSubgroup))])

%% Quotient Group/商群
if ~isempty(NormalSubgroup)
    FactorGroup = Square_group/NormalSubgroup{3}; % Create quotient group/创建商群
    disp('Quotient group structure:')
    disp(FactorGroup)
end

%% Group Generators/群生成元
[Generators] = Square_group.generator(); % Find minimal generators/寻找最小生成元
disp(['Group rank: ' num2str(length(Generators{1}))])