%% Initialize the Hopping Parameters
% We start by defining symbolic variables and setting up the hopping parameters for a graphene lattice.
% Initialize symbolic variable for hopping parameter
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
%% Print Out the Hopping Information
% Print out current Hamiltonian configuration
Graphene_TB.printout;
%% Convert and Display the List Representation
% Convert to list representation
Graphene_TB_list = Graphene_TB.rewrite();
Graphene_TB_list.list();
%% Convert Back to Matrix Representation
% Convert back to matrix form
Graphene_TB_list2mat = Graphene_TB_list.rewind();
Graphene_TB_list2mat.printout;
%% Set Up the Lattice
% We define the lattice structure, either by directly using a POSCAR file or defining vectors manually.
% Input orbital structure from a POSCAR file
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

% Further Manipulations
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
%% Overload test
