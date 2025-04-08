function [Rnn,nn_store_smart,Atom_store_smart,Rnn_map] = nn_information(H_hr,silence)
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
