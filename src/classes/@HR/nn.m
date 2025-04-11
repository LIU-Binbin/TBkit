function H_hr = nn(H_hr,search_range,Accuracy,Rlength_cut,options)
% NN Generate nearest neighbor information for HR object
%
%   H_HR = NN(H_HR,SEARCH_RANGE,ACCURACY,RLENGTH_CUT,OPTIONS) generates
%   nearest neighbor information for the HR object
%
%   Inputs:
%       H_hr - HR object to process
%       search_range - Search range vector [default: [0 0 0]]
%       Accuracy - Rounding accuracy [default: 1e-4]
%       Rlength_cut - Cutoff distance [default: 5]
%       options.MAX_NN_LENGTH - Maximum neighbor count [default: 1e7]
%       options.onsite - Flag for onsite terms [default: false]
%   Output:
%       H_hr - Updated HR object with neighbor information
%
%   Notes:
%       - Inherits from TBkit class implementation
%       - Uses named arguments for optional parameters
arguments
H_hr HR;
search_range double = [0 0 0];
Accuracy double = 1e-4;
Rlength_cut double = 5;
options.MAX_NN_LENGTH =  10000000;
options.onsite = false;
end
optionsCell = namedargs2cell(options);
H_hr = nn@TBkit(H_hr,search_range,Accuracy,Rlength_cut,optionsCell{:});
end
