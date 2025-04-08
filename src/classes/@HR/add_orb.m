function H_hr = add_orb(H_hr,hop_struct,orbOne,QuantumOne,elementOne)
nstruct = length( hop_struct);
if nargin < 3
orbOne = repmat([0,0,0],[nstruct,1]);
end
if nargin < 4
QuantumOne = repmat([1,0,0,1],[nstruct,1]);
end
if nargin < 5
elementOne = ones(nstruct,1);
end
WANNUM = H_hr.WAN_NUM;
H_hr = H_hr.expand_empty_one(orbOne,QuantumOne,elementOne);
for j = 1:nstruct
for i = 1:length(hop_struct(j).hop)
H_hr = H_hr.set_hop(hop_struct(j).hop(i),hop_struct(j).hi(i),WANNUM+j,hop_struct(j).vector,'set');
H_hr = H_hr.set_hop(conj(hop_struct(j).hop(i)),WANNUM+j,hop_struct(j).hi(i),-hop_struct(j).vector,'set');
end
end
end
