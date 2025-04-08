function H_hr = add_orb(H_hr, hop_struct, orbOne, QuantumOne, elementOne)
%ADD_ORB Expands Hamiltonian with new orbitals and hopping terms
%   This function extends existing tight-binding Hamiltonian by adding new
%   orbitals and their corresponding hopping terms while preserving Hermiticity
%
%   Inputs:
%       H_hr       - Initial Hamiltonian object with WAN_NUM property
%       hop_struct - Array of structs containing hopping parameters:
%                    .hop   : Hopping amplitudes (complex numbers)
%                    .hi    : Target orbital indices
%                    .vector: Lattice vectors for hopping [3x1]
%       orbOne     - (Optional) Orbital positions [nstruct x 3], default zeros
%       QuantumOne - (Optional) Quantum numbers [nstruct x 4], default [1,0,0,1]
%       elementOne - (Optional) Element type indices [nstruct x 1], default 1
%
%   Output:
%       H_hr       - Updated Hamiltonian with new orbitals and hoppings
%
%   Features:
%       - Automatically adds Hermitian conjugate terms
%       - Supports multiple orbital configurations
%       - Default parameters for rapid prototyping

    nstruct = length(hop_struct);
    
    % Set default orbital positions if not provided
    if nargin < 3
        orbOne = repmat([0,0,0], [nstruct,1]);
    end
    
    % Set default quantum numbers if not provided
    if nargin < 4
        QuantumOne = repmat([1,0,0,1], [nstruct,1]);
    end
    
    % Set default element types if not provided
    if nargin < 5
        elementOne = ones(nstruct,1);
    end
    
    % Get current number of orbitals
    WANNUM = H_hr.WAN_NUM;
    
    % Expand Hamiltonian basis with new orbitals
    H_hr = H_hr.expand_empty_one(orbOne, QuantumOne, elementOne);
    
    % Add hopping terms for each structure configuration
    for j = 1:nstruct
        for i = 1:length(hop_struct(j).hop)
            % Add forward hopping term
            H_hr = H_hr.set_hop(hop_struct(j).hop(i), ...
                hop_struct(j).hi(i), WANNUM+j, hop_struct(j).vector, 'set');
            
            % Add Hermitian conjugate term (reverse hopping)
            H_hr = H_hr.set_hop(conj(hop_struct(j).hop(i)), ...
                WANNUM+j, hop_struct(j).hi(i), -hop_struct(j).vector, 'set');
        end
    end
end