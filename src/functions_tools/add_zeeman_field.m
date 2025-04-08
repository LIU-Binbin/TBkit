function Ham_obj_B = add_zeeman_field(Ham_obj, B, opts)
% ADD_ZEEMAN_FIELD Adds a Zeeman field to the Hamiltonian.
% The function modifies the Hamiltonian object by applying the Zeeman effect
% due to an applied magnetic field B.
%
% The Zeeman term is given by:
%   H_Zeeman = -mu_B * g_factor * B * S,
% where:
%   mu_B is the Bohr magneton,
%   g_factor is the g-factor of the material,
%   B is the magnetic field (Bx, By, Bz),
%   S is the spin operator.
%
% Input:
%   Ham_obj - An HR object containing the Hamiltonian data (HnumL, Nb, etc.).
%   B - Magnetic field vector (default: [1 0 0] Tesla).
%   opts - A structure of options:
%     - opts.S: Spin quantum number (default: 1/2).
%     - opts.L: Orbital angular momentum (default: 0).
%     - opts.J: Total angular momentum (default: opts.L + opts.S).
%     - opts.g_factor: g-factor (default: 2).
%     - opts.spin_mode: Spin ordering ('new' or 'old', default: 'new').
%     - opts.natural_unit: Use natural units (default: false).
%
% Output:
%   Ham_obj_B - The modified Hamiltonian object with the Zeeman term added.

    % Validate input arguments
    arguments
        Ham_obj
        B double = [1 0 0]  % Tesla, Magnetic field components [Bx, By, Bz]
        opts.S double = 1/2  % Spin quantum number (default 1/2)
        opts.L double = 0    % Orbital angular momentum (default 0)
        opts.J double = 0    % Total angular momentum (default opts.L + opts.S)
        opts.g_factor double = 2  % g-factor (default 2)
        opts.spin_mode {mustBeMember(opts.spin_mode, {'new', 'old'})} = 'new'  % Spin ordering
        opts.natural_unit = false  % Use natural units for mu_B (default: false)
    end

    % Ensure that the Hamiltonian object has numerical data (HnumL)
    if isempty(Ham_obj.HnumL)
        error("Please input a valid numerical HR object.");
    end

    % Pauli matrices for spin operators
    sigma_x = [0 1; 1 0];
    sigma_y = [0 -1i; 1i 0];
    sigma_z = [1 0; 0 -1];

    % Magnetic moment (mu_B) in eV/Tesla or use natural units if specified
    if opts.natural_unit
        mu_B = 1;  % Use natural units (no conversion)
    else
        mu_B = 5.7884e-5;  % eV/Tesla (physical constant)
    end

    % Compute the dot product of the magnetic field and the spin operators
    s_dot_B = mu_B * (B(1) * sigma_x + B(2) * sigma_y + B(3) * sigma_z);

    % Set the total angular momentum J if not provided
    Nproj = Ham_obj.Nbands / 2;  % Number of projections (assuming spin-1/2)
    if opts.J == 0
        opts.J = opts.L + opts.S;  % Default J is L + S
    end

    % Identity matrix for the projection part of the Hamiltonian
    I_N_proj = diag(ones(1, Nproj) * opts.J * opts.g_factor);

    % Apply the spin mode selection (new or old ordering)
    if opts.spin_mode == "new"
        disp("Spin indexes of projs = 1-up, 1-dn, 2-up, 2-dn, ...");
        disp("This is used in Wannier90 v2.x and newer.");
        % Construct Hamiltonian in the new spin ordering: kron(I_N_proj, s_dot_B)
        H_B = kron(I_N_proj, s_dot_B);
    elseif opts.spin_mode == "old"
        disp("Spin indexes of projs = 1-up, 2-up, ... 1-dn, 2-dn, ...");
        disp("This is used in Wannier90 v1.2 and Wannier Tools.");
        % Construct Hamiltonian in the old spin ordering: kron(s_dot_B, I_N_proj)
        H_B = kron(s_dot_B, I_N_proj);
    end

    % Add the Zeeman term to the original Hamiltonian
    Ham_obj_B = Ham_obj + H_B;  % Final Hamiltonian with Zeeman field

end

% function Ham_obj_B = add_zeeman_field(Ham_obj,B,opts)
% % H_B = g*mu_B/hbar*s*B = mu_B*sigma*B (default: g=2, s=hbar/2*sigma)
% arguments
%     Ham_obj
%     B double = [1 0 0] % Tesla, [Bx By Bz]
%     opts.S double = 1/2
%     opts.L double = 0
%     opts.J double = 0
%     opts.g_factor double = 2
%     opts.spin_mode {mustBeMember(opts.spin_mode,{'new','old'})} = 'new'
%     opts.natural_unit = false
% end
% %%
% if isempty(Ham_obj.HnumL)
%     error("Please input a numerical HR obj")
% end
% %%
% sigma_x = [0 1;1 0]; sigma_y = [0 -1i;1i 0]; sigma_z = [1 0;0 -1];
% if opts.natural_unit
%     mu_B = 1;
% else
%     mu_B = 5.7884e-5; % eV/Tesla
% end
% s_dot_B = mu_B.*(B(1)*sigma_x + B(2)*sigma_y + B(3)*sigma_z);
% %%
% Nproj = Ham_obj.Nbands/2;
% if opts.J == 0
%     opts.J = opts.L + opts.S;
% end
% I_N_proj = diag(ones(1,Nproj) .* opts.J .* opts.g_factor);
% 
% if opts.spin_mode == "new"
%     disp("spin indexes of projs = 1-up, 1-dn, 2-up, 2-dn, ...")
%     disp("This is used in Wannier90 v2.x and newer")
%     H_B = kron(I_N_proj, s_dot_B);
% elseif opts.spin_mode == "old"
%     disp("spin indexes of projs = 1-up, 2-up, ... 1-dn, 2-dn, ...")
%     disp("This is used in Wannier90 v1.2 and Wannier Tools")
%     H_B = kron(s_dot_B, I_N_proj);
% end
% Ham_obj_B = Ham_obj + H_B;
% end