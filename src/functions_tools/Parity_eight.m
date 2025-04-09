function parityList = Parity_eight(occupancyData, labelIndices)
% PARITY_EIGHT Calculate parity values for multiple label configurations
%   This function computes parity values for a series of label files using
%   the provided occupancy data. Supports custom label index ranges.
%
% Inputs:
%   occupancyData  - Input occupancy data matrix/array
%   labelIndices   - (Optional) Array of label indices to process. 
%                    Default: 1:8
%
% Output:
%   parityList     - Array of computed parity values corresponding to
%                    each label index
%
% Example:
%   parityResults = Parity_eight(occupancyMatrix, [1 3 5 7])

% Handle optional input
if nargin < 2
    labelIndices = 1:8;
end

% Validate input format
validateattributes(labelIndices, {'numeric'}, {'vector', 'integer', 'positive'}, ...
                   'Parity_eight', 'labelIndices', 2);

% Preallocate output array for performance
% parityList = zeros(size(labelIndices));

% Process each label index
for idx = 1:numel(labelIndices)
    currentIndex = labelIndices(idx);
    
    % Generate formatted label filename
    labelFile = sprintf('label_file_%d', currentIndex);
    
    % Calculate parity and store result
    parityList(idx) = Parity(labelFile, occupancyData);
end
end