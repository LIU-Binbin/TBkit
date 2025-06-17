function [nn_sparse_temp, Rnn_list] = nn_sparse_gen(orb1, orb2, Rm, search_rangex, search_rangey, search_rangez, Accuracy, Rlength_cut, onsite)
%NN_SPARSE_GEN Generate sparse neighbor list for tight-binding models
%   This function calculates valid hopping terms between two orbitals within specified
%   search ranges and distance cutoff. Supports both numerical and symbolic computations.
%
%   Inputs:
%     orb1, orb2      : Fractional coordinates of two orbitals in unit cell (1x3 vectors)
%     Rm              : Lattice vectors matrix (3x3) for converting fractional to Cartesian
%     search_range[xyz]: Search range in each lattice direction (integers)
%     Accuracy        : Decimal precision for distance rounding (-1 for symbolic mode)
%     Rlength_cut     : Maximum hopping distance cutoff
%     onsite          : Flag to include same-site terms (default = false)
%
%   Outputs:
%     nn_sparse_temp  : Sparse matrix containing hopping info [i,j,dx,dy,dz,Rx,Ry,Rz,dist,hop]
%     Rnn_list        : List of neighbor distances

% Initialize parameters
Rc1 = orb1;
Rc2 = orb2;
R_fractional_diff = -(Rc1 - Rc2);  % Relative position vector
count = 1;

% Preallocate arrays
reducible_num = (2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
Rnn_list = zeros(reducible_num, 1);

% Initialize output matrix based on computation mode
symmode = (Accuracy == -1);
if symmode
    nn_sparse_temp = sym(zeros(reducible_num, 10));
else
    nn_sparse_temp = zeros(reducible_num, 10);
end
AccuracyNum = round(log(Accuracy)/log(10));
% Main search loop over neighboring cells
for Rf_a1 = -search_rangex:search_rangex
    for Rf_a2 = -search_rangey:search_rangey
        for Rf_a3 = -search_rangez:search_rangez
            R_vector = [Rf_a1, Rf_a2, Rf_a3];
            Rij_cart = (R_vector + R_fractional_diff) * Rm;
            Rlength = norm(Rij_cart);

            % Apply numerical rounding or symbolic comparison
            if ~symmode
                Rlength = round(Rlength, -AccuracyNum);
                valid = ((0 < Rlength) || (Rlength == 0 && onsite)) && (Rlength < Rlength_cut);
            else
                valid = (sym(0) < Rlength || (sym(0) == Rlength && onsite)) && (Rlength < sym(Rlength_cut));
            end

            % Store valid neighbor information
            if valid
                nn_sparse_temp(count, 6:8) = R_vector;
                nn_sparse_temp(count, 3:5) = Rij_cart;
                nn_sparse_temp(count, 9) = Rlength;
                Rnn_list(count) = Rlength;
                count = count + 1;
            end
        end
    end
end

% Trim unused preallocated space
if count <= reducible_num
    Rnn_list(count:end) = [];
    nn_sparse_temp(count:end, :) = [];
end
end

function y = roundn(x, n)
%ROUNDN Round to nearest 10^n (compatibility function)
y = round(x*(10^-n))*(10^n);
end