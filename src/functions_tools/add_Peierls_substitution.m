function add_Peierls_substitution(Ham_obj)
% ADD_PEIERLS_SUBSTITUTION Modifies the Hamiltonian matrix with the Peierls substitution.
% This function applies a magnetic field through the Peierls phase to the Hamiltonian.
% The Peierls substitution is implemented via the A (vector potential) field.
%
% Input:
%   Ham_obj - An object of class HR, which contains the Hamiltonian matrix and necessary parameters.
%             The object should have the following fields:
%             - num: Numerical flag (non-zero if valid)
%             - vectorL: Vector array for lattice points
%             - orbL: Orbital positions in the lattice
%             - Rm: Lattice scaling factor
%             - WAN_NUM: Number of Wannier functions
%             - NRPTS: Number of grid points
%             - HnumL: Hamiltonian matrix (Lattice version)
%
% Notes:
% - A magnetic field is applied along the x-axis (B = 1T along x), 
%   and the phase factor is added to the Hamiltonian.
% - The value of phi0 is the magnetic flux quantum.

    % Validate input argument Ham_obj
    arguments
        Ham_obj HR  % HR class object with necessary fields
    end

    % Check if the Ham_obj contains a valid numerical object
    if ~Ham_obj.num
        error("Please input a valid numerical HR object");
    end

    % Constants
    A = [1 0 0];  % Magnetic field vector (B = 1T along x-direction)
    phi0 = 4.13567e-15;  % Magnetic flux quantum in Wb·T·m²
    phi0 = phi0 * 1e20;  % Convert to Wb·T·Å² (angstrom^2)

    % Calculate the transformed lattice vectors and orbital positions
    vectorL_cart = Ham_obj.vectorL * Ham_obj.Rm;  % Lattice vectors in Cartesian coordinates
    orbL_cart = Ham_obj.orbL * Ham_obj.Rm;        % Orbital positions in Cartesian coordinates

    % Loop over all Wannier functions (i, j) and grid points (k)
    for i = 1:Ham_obj.WAN_NUM
        for j = 1:Ham_obj.WAN_NUM
            for k = 1:Ham_obj.NRPTS
                % Calculate the mid-point between the two orbitals and the k-th lattice vector
                xyz_mid = ( orbL_cart(i,:) + orbL_cart(j,:) + vectorL_cart(k,:) ) / 2;  % Å
                
                % Compute the vector potential A field at the midpoint
                A_field = A * xyz_mid';  % Vector potential: (y * B, 0, 0) in Cartesian coordinates
                
                % Calculate the Peierls phase factor based on the displacement between orbitals
                A_dl = A_field * ( orbL_cart(j,1) + vectorL_cart(k,1) - orbL_cart(i,1) );  % Displacement in x-direction

                % Apply the Peierls phase to the Hamiltonian matrix element
                Ham_obj.HnumL(i,j,k) = Ham_obj.HnumL(i,j,k) * exp(1i * 2*pi / phi0 * A_dl);
            end
        end
    end

end


% function add_Peierls_substitution(Ham_obj)
% % H_B = g*mu_B/hbar*s*B = mu_B*sigma*B (default: g=2, s=hbar/2*sigma)
% arguments
%     Ham_obj HR
% end
% %%
% if ~Ham_obj.num
%     error("Please input a numerical HR obj")
% end
% %%
% vectorL_cart = Ham_obj.vectorL * Ham_obj.Rm;
% orbL_cart = Ham_obj.orbL * Ham_obj.Rm;
% A = [1 0 0]; % B=1T
% phi0 = 4.13567e-15; % Wb, T*m^2
% phi0 = phi0 * 1e20; % Wb, T*Ang^2
% %%
% for i = 1:Ham_obj.WAN_NUM
%     for j = 1:Ham_obj.WAN_NUM
%         for k = 1:Ham_obj.NRPTS
%             xyz_mid = ( orbL_cart(i,:) + orbL_cart(j,:) + vectorL_cart(k,:) )./2; % Ang
% 
%             A_field = A * xyz_mid'; % (y*B, 0 , 0)
%             A_dl = A_field * ( orbL_cart(j,1) + vectorL_cart(k,1) - orbL_cart(i,1) ); % x_j - x_i
% 
%             Ham_obj.HnumL(i,j,k) = Ham_obj.HnumL(i,j,k) * exp(1i * 2*pi/phi0 * A_dl);
%         end
%     end
% end
% 
% end