function [Rnn,nn_store_smart,Atom_store_smart,Rnn_map] = nn_information(H_hr,silence)
% NN_INFORMATION Retrieve nearest neighbor information from HR object
%
%   [RNN,NN_STORE_SMART,ATOM_STORE_SMART,RNN_MAP] = NN_INFORMATION(H_HR,SILENCE)
%   extracts stored nearest neighbor information
%
%   Inputs:
%       H_hr - HR object containing NN information
%       silence - Boolean to suppress output [default: true]
%   Outputs:
%       Rnn - Array of neighbor distances
%       nn_store_smart - Smart neighbor storage structure
%       Atom_store_smart - Smart atom storage structure
%       Rnn_map - Container mapping distances to indices
%
%   Notes:
%       - Provides both smart and regular NN storage access
%       - Optionally displays neighbor distances when not silent
if nargin < 2
silence = true;
end
Rnn = H_hr.Rnn;
if ~isempty(H_hr.nn_store_smart)
nn_store_smart = H_hr.nn_store_smart;
else
nn_store_smart = H_hr.nn_store;
end
try
Atom_store_smart = H_hr.Atom_store_smart;
catch
Atom_store_smart = [];
end
Rnn_map = H_hr.Rnn_map ;
if ~silence
for i = 1:length(Rnn)
fprintf('the %3d th Rnn vector: %7.4f Angstrom\n',i,Rnn(i));
end
end
end
