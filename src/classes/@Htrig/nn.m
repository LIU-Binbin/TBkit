% NN Compute nearest neighbor information for a Htrig object.
%
% SYNTAX:
%   H_htrig = nn(H_htrig, search_range, Accuracy, Rlength_cut)
%
% DESCRIPTION:
%   This method overloads the nn function from the TBkit superclass for the Htrig
%   object. It computes the nearest neighbor information based on the provided search range,
%   accuracy tolerance, and a cutoff value for the distance R. If these parameters are not
%   supplied, default values are used:
%       - search_range defaults to [0, 0, 0].
%       - Accuracy defaults to 4.
%       - Rlength_cut defaults to 5.
%
%   After computing the nearest neighbors using the superclass method, the Htrig object's
%   Type is set to 'exp' to reflect the updated representation.
%
% INPUTS:
%   H_htrig     - A Htrig object.
%   search_range- A vector specifying the range to search for nearest neighbors (default: [0,0,0]).
%   Accuracy    - A numeric value specifying the accuracy tolerance for the neighbor search (default: 4).
%   Rlength_cut - A numeric value specifying the cutoff length for R (default: 5).
%
% OUTPUT:
%   H_htrig     - The updated Htrig object with computed nearest neighbor information.
%
% EXAMPLE:
%   % Compute nearest neighbors with default parameters:
%   H_htrig = nn(H_htrig, [0,0,0], 4, 5);
%
function H_htrig = nn(H_htrig, search_range, Accuracy, Rlength_cut)
    if nargin < 4
        Rlength_cut = 5;
    end
    if nargin < 3
        Accuracy = 4;
    end
    if nargin < 2
        search_range = [0 0 0];
    end
    H_htrig = nn@TBkit(H_htrig, search_range, Accuracy, Rlength_cut);
    H_htrig.Type = 'exp';
end

