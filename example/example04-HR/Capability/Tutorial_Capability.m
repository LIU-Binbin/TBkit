%% Tutorial: Capability of HR class
%% Initialize the Hopping Parameters
% We start by defining symbolic variables and setting up the hopping parameters 
% for a graphene lattice. Initialize symbolic variable for hopping parameter

syms t real;

% Create a tight-binding Hamiltonian representation for graphene
Graphene_TB = HR(2);

% Define nearest-neighbor vectors for hopping
A_vectorL = [0,0,0; 1,0,0; 0,-1,0];
B_vectorL = [0,0,0; -1,0,0; 0,1,0];

% Set up hopping terms between lattice sites
Graphene_TB = Graphene_TB.set_hop(t, 1, 2, A_vectorL(1,:), 'sym');
Graphene_TB = Graphene_TB.set_hop(t, 1, 2, A_vectorL(2,:), 'sym');
Graphene_TB = Graphene_TB.set_hop(t, 1, 2, A_vectorL(3,:), 'sym');

% Set multivector hopping for additional connections
Graphene_TB = Graphene_TB.set_hop(t, 2, 1, B_vectorL, 'sym');
% Print Out the Hopping Information
% Print out current Hamiltonian configuration

Graphene_TB.printout;
% Convert and Display the List Representation
% Convert to list representation

Graphene_TB_list = Graphene_TB.rewrite();
Graphene_TB_list.list();
% Convert Back to Matrix Representation
% Convert back to matrix form

Graphene_TB_list2mat = Graphene_TB_list.rewind();
Graphene_TB_list2mat.printout;
%% Set Up the Lattice
% We define the lattice structure, either by directly using a POSCAR file or 
% defining vectors manually. Input orbital structure from a POSCAR file

Graphene_TB = Graphene_TB.input_orb_struct('POSCAR', 'tbsk');

% Alternatively, use the overload operator '<' for importing
Graphene_TB = Graphene_TB < 'POSCAR';

% Direct input of lattice vectors
Rm = [1,0,0; -0.5, sqrt(3)/2, 0; 0,0,1];
Graphene_TB_list = Graphene_TB_list.input_Rm(Rm);

% Define orbital positions in the lattice
orbL = [2/3, 1/3, 0; 1/3, 2/3, 0];
Graphene_TB_list.orbL = orbL;

% Display the hopping structure based on the setup
Graphene_TB_list.show('HOPPING', 'scale', 2.4560000896, 'atomscale', 1, 'TwoD', true);
%% Convert Hamiltonian to Bloch representation and simplify
% This convention simplifies the Hamiltonian without removing positional information.

H_Graphene_TB_convention_I = simplify(Graphene_TB_list.sym());
%  Omit the Ionic Positions
H_Graphene_TB_convention_II = simplify(Graphene_TB_list.sym("simple",1))
%  Use Relative Lattice Vectors
H_Graphene_TB_tmp = simplify(Graphene_TB_list.sym("simple",1,'cartesian',0))

%% Further Manipulations
% Example use of MATLAB's 'latex' function to convert symbolic expressions
latex_representation = latex(H_Graphene_TB_convention_I)
% Use 'rewrite' function for additional forms, such as trigonometrics or others
rewritten_expression = rewrite(H_Graphene_TB_convention_I, 'sincos')
%% Taylor Expansion
% Define the fractional coordinates for the k-point K

K_frac = [1/3, 1/3, 0];
% Convert HR object to HK object with a Taylor expansion around K_frac
H_Graphene_K = Graphene_TB_list.HR2HK(K_frac, "sym", 1, "Order", 1);
% Declare symbolic variables for k_x and k_y
syms k_x k_y;
% Simplify the Hamiltonian and collect terms
H_Graphene_K_simplified = collect(simplify(H_Graphene_K.sym, "Steps", 1000), k_x - 1i*k_y);
% Define the fractional coordinates for the k-point K'
K_frac_prime = -K_frac;
% Convert the HR object to HK object for the  k-point  K'
H_Graphene_K__prime = Graphene_TB_list.HR2HK(K_frac_prime, "sym", 1, 'Order', 2);
% Simplify the Hamiltonian at K' and collect terms
collect(simplify(H_Graphene_K__prime.sym,"Steps",1000),k_x - 1i*k_y)
% Finally, you can obtain the conjugate of the Hamiltonian at the K' is same as K, which is a common practice in solid-state physics.
% H_K = conj(H_Graphene_K__prime);
%% Tutorial: Operator Overloading in HR Classes
% This tutorial demonstrates operator overloading implementation for a custom 
% graphene tight-binding (TB) class.
%% 1. Overloading Basic Arithmetic Operators
% The '+' operator is overloaded to handle both homogeneous and heterogeneous 
% object types.

% 1.1 Adding two Graphene_TB objects
obj_sum = Graphene_TB + Graphene_TB;
obj_sum.printout();  % Display resultant object properties

% 1.2 Mixed-type addition
half_list = Graphene_TB_list.reseq(':',1:3);  % Create subset
combined = half_list + Graphene_TB_list;      % List-to-list addition
combined.list();                              % Display list contents
%% 2. Compound Operations and Type Handling
% 2.1 Cross-type operations

mixed_result = Graphene_TB + Graphene_TB_list;
mixed_result.printout();

% 2.2 Chained operations
final_obj = Graphene_TB + Graphene_TB_list - Graphene_TB;
final_obj.printout();
%% 3. Element-wise and Matrix Operations
% 3.1 Element-wise multiplication (.*)

A = diag([2,1]) + diag([-1],-1);         % Create test matrix
elementwise_result = A .* Graphene_TB;   % Requires matrix-type support
elementwise_result.printout();

% 3.2 Scalar multiplication
scaled_obj = 2 .* Graphene_TB;           % Scalar expansion
scaled_obj.printout();

% 3.3 Matrix multiplication (*)
unitary_op = [0 1i; -1i 0];              % Pauli-Y matrix
transformed = unitary_op * Graphene_TB * unitary_op';
transformed.printout();
%% 4. Power Operators and Complex Operations
% 4.1 Element-wise power

squared_elementwise = Graphene_TB.^2;
squared_elementwise.printout();

% 4.2 Matrix power
matrix_power = Graphene_TB^2;            % Requires square matrix
matrix_power.printout();
%% 5. Matrix Transformations and Symbolic Verification
% 5.1 Create symbolic TB Hamiltonian

syms k_x k_y real;
H_TB = Graphene_TB_list.rewind() + 1i*Graphene_TB_list.rewind();
H_sym = H_TB.sym('simple',1);           % Symbolic simplification
H_TB.printout();

% 5.2 Transposition operations
% Regular transpose (.')
H_transposed = H_TB.';
H_transposed_sym = H_transposed.sym('simple',1);
disp('Transposition difference:');
simplify(H_sym.' - H_transposed_sym)    % Should return zero

% Hermitian transpose (')
H_hermitian = H_TB';
H_herm_sym = H_hermitian.sym('simple',1);
disp('Conjugate transpose difference:');
simplify(H_sym' - H_herm_sym)           % Verify correctness

% 5.3 Complex conjugation
H_conj = conj(H_TB);
H_conj_sym = H_conj.sym('simple',1);
disp('Conjugation difference:');
simplify(conj(H_sym) - H_conj_sym)      % Validate conjugation
%% Tutorial: Kronecker Product Operations in Spinful Graphene Models
%% 1. Constructing Spinful Hamiltonians
% Use MATLAB's `kron` to create spinful configurations:

% 1.1 Block-diagonal spin structure (uu/dd basis)
Graphene_spinful = [Graphene_TB,Graphene_TB];
Graphene_spinful = kron(eye(2), Graphene_TB);  % 2×2 identity ⊗ TB model
disp('Block-diagonal spin structure:');
Graphene_spinful.sym("simple", 1, 'cartesian', 0)  % Symbolic simplification
Graphene_spinful.orbL

% 1.2 Orbital-interleaved structure (ud/ud basis)
Graphene_interleaved = kron(Graphene_TB, eye(2));  % Reverse Kronecker order
disp('Interleaved spin-orbit structure:');
Graphene_interleaved.sym("simple", 1, 'cartesian', 0)
Graphene_interleaved.orbL
%% 2. Index Resequencing for Spin Configuration Control
% 2.1 Custom orbital ordering

Graphene_spinful = Graphene_interleaved.reseq([1,3,2,4]);  % Reorder indices
disp('Custom uudd ordering:');
Graphene_spinful.sym("simple", 1, 'cartesian', 0)

% 2.2 Construct spinless subsystem
Graphene_spinless = Graphene_interleaved.reseq([1,3]);  % Select subset
disp('Spinless reduction:');
Graphene_spinless.sym("simple", 1, 'cartesian', 0)
%% 3. Symbolic Workflow Verification

syms k_x k_y real;
disp('Symbolic equivalence check:');
% 3.1 Verify Kronecker properties
H_block = kron(eye(2), Graphene_TB.sym());
H_inter = kron(Graphene_TB.sym(), eye(2));
disp('Block vs. interleaved trace difference:');
simplify(trace(H_block) - 2*trace(Graphene_TB.sym()))  % Should be 0
disp('Interleaved vs. original dimension ratio:');
size(H_inter,1)/size(Graphene_TB.sym(),1)  % Should be 2
% 3.2 Validate index resequencing
original_dim = size(Graphene_spinful.sym(), 1);
reduced_dim = size(Graphene_spinless.sym(), 1);
disp(['Dimension reduction factor: ' num2str(original_dim/reduced_dim)]);
% - `kron` order determines spin-orbital coupling structure:
%   kron(eye(n), A) -> Block diagonal     (n copies)
%   kron(A, eye(n)) -> Element-wise scale (interleaved)
%
% - `reseq` parameters: [new_order] specifies orbital indices permutation
% - `cartesian` 0 flag: Use fractional coordinates for k-space representation
% - Symbolic workflow requires consistent variable declarations (k_x/k_y)
%% Tutorial: Band Structure Calculation for Graphene
% This tutorial demonstrates numerical band structure computation for monolayer graphene
% using HR class implementation.
t = 1;  % Nearest-neighbor hopping energy [eV]
%%  Brillouin Zone Path Generation
% Load k-path from file (alternative method)
%  Graphene_TB = Graphene_TB < 'KPOINTS';
% Graphene_TB = Graphene_TB.kpathgen3D('KPOINTS');

% 2.2 Manual k-path specification
high_sym_points = [  
    0.0000000000   0.0000000000   0.0000000000     ;  % Γ point
    0.5000000000   0.0000000000   0.0000000000     ;  % M point
    0.5000000000   0.0000000000   0.0000000000     ;  % M point
    0.3333333333   0.3333333333   0.0000000000     ;  % K point
    0.3333333333   0.3333333333   0.0000000000     ;  % K point
    0.0000000000   0.0000000000   0.0000000000     ;  % Γ point (return)    
];

segment_nodes = [60, 50, 30];  % Points per path segment [Γ-M, M-K, K-Γ]

% Generate k-path coordinates
[Graphene_TB.klist_cart, Graphene_TB.klist_frac, Graphene_TB.klist_l, ...
 Graphene_TB.kpoints_l, Graphene_TB.kpoints_frac] = ...
    TBkit.kpathgen(high_sym_points, segment_nodes, Graphene_TB.Gk);

% Set high-symmetry point labels
Graphene_TB.kpoints_name = ["Γ", "M", "K", "Γ"];

% Convert symbolic parameters to numerical values
Graphene_TB_numerical = Graphene_TB.Subsall();  % Substitute symbolic variables

% Generate eigenvalue spectrum
EIGENCAR = Graphene_TB_numerical.EIGENCAR_gen();  % Compute band structure


%  Object-oriented plotting (uses internal k-path data)
Graphene_TB_numerical.bandplot(EIGENCAR, [-3, 3],...
    'Title', 'Graphene Band Structure (t = 1 eV)',...
    'YLabel', 'Energy [eV]');

%  Direct function call equivalent
bandplot(EIGENCAR, [-3, 3],...
    Graphene_TB_numerical.klist_l,...
    Graphene_TB_numerical.kpoints_l,...
    Graphene_TB_numerical.kpoints_name,...
    'Title', 'Graphene Band Structure (t = 1 eV)',...
    'YLabel', 'Energy [eV]');
%
% - `Subsall` replaces symbolic variables with numerical values
% - `EIGENCAR_gen` diagonalizes Hamiltonian at each k-point
% - `bandplot` handles both object-contained and explicit k-path data
%% reshape 
[~,Ax] = Figs(4,2);
%
Graphene_TB.show('TwoD',1,'ax',Ax(1),'cmap',@hsv);title(Ax(1),'Primitive')
Graphene_TB.bandplot('ax',Ax(5),'title','Primitive');
%
Graphene_TB_221 = Graphene_TB.supercell_hr(diag([2,2,1]));
Graphene_TB_221.show('TwoD',1,'ax',Ax(2),'cmap',@hsv);title(Ax(2),'Supercell 221')
Graphene_TB_221.bandplot('ax',Ax(6),'title','Supercell 221');
%
Ns = [1 0 0;1 2 0;0 0 1];
Graphene_TB_rec = Graphene_TB.supercell_hr(Ns); 
Graphene_TB_rec = Graphene_TB_rec < 'KPOINTS_rec';
Graphene_TB_rec.show('TwoD',1,'ax',Ax(3),'cmap',@hsv);
title(Ax(3),'Supercell rec')
Graphene_TB_rec.bandplot('ax',Ax(7),'title','Supercell rec');
%
Ns = [1 -1 0;1 2 0;0 0 1];
Graphene_TB_arm = Graphene_TB.supercell_hr(Ns);
Graphene_TB_arm.show('TwoD',1,'ax',Ax(4),'cmap',@hsv);title(Ax(4),'Supercell arm')
Graphene_TB_arm.bandplot('ax',Ax(8),'title','Supercell arm');
