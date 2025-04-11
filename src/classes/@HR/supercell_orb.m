function [sc_orb,sc_vec,sc_elementL,sc_quantumL] = supercell_orb(H_hr,Ns,Accuracy)
% SUPERCELL_ORB Generate supercell orbital positions
%
%   [SC_ORB,SC_VEC,SC_ELEMENTL,SC_QUANTUML] = SUPERCELL_ORB(H_HR,NS,ACCURACY)
%   calculates supercell orbital positions
%
%   Inputs:
%       H_hr - Original HR object
%       Ns - Supercell transformation matrix [default: eye(Dim)]
%       Accuracy - Rounding accuracy [default: 1e-6]
%   Outputs:
%       sc_orb - Supercell orbital positions
%       sc_vec - Supercell vectors
%       sc_elementL - Element list
%       sc_quantumL - Quantum numbers
%
%   Notes:
%       - Validates supercell matrix
%       - Handles fractional coordinates
%       - Maintains periodic boundary conditions
arguments
H_hr HR;
Ns double = eye(H_hr.Dim);
Accuracy double = 1e-6;
end
if ~(Ns == round(Ns))
error("sc_red_lat array elements must be integers");
end
if det(Ns) < Accuracy
error("Super-cell lattice vectors length/area/volume too close to zero, or zero.");
end
if det(Ns)<0.0
error("Super-cell lattice vectors need to form right handed system.");
end
nAccuracy = round(log10(Accuracy));
orb_init = H_hr.orbL;
DIM = H_hr.Dim;
max_R=max(abs(Ns))*3;
for d = 1:DIM
dL{d} = (-max_R(d):max_R(d)).';
end
sc_cands = fold(@park.SemiProductVector,dL);
sc_vec=([]);
eps_shift=sqrt(2.0)*1.0E-8;
for ivec = 1:length(sc_cands)
tmp_red = TBkit.to_red_sc(sc_cands(ivec,:),Ns);
inside = 1;
for it = 1:length(tmp_red)
t = tmp_red(it);
if t <= -1.0*eps_shift || t>1.0-eps_shift
inside=0;
end
end
if inside == 1
sc_vec=[sc_vec;sc_cands(ivec,:)];
end
end
[num_sc,~] = size(sc_vec);
if round(round(abs(det(Ns)))) ~= num_sc
error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
end
WANNUM = H_hr.WAN_NUM;
sc_orb=zeros(num_sc*WANNUM,H_hr.Dim);
sc_elementL = repmat(H_hr.elementL,[num_sc,1]);
sc_quantumL = repmat(H_hr.quantumL,[num_sc,1]);
count_tmp = 1;
for icur_sc_vec = 1:num_sc
cur_sc_vec = sc_vec(icur_sc_vec,:);
for iorb = 1:WANNUM
sc_orb(count_tmp,:) = roundn((orb_init(iorb,:)+cur_sc_vec)/Ns,nAccuracy);
count_tmp =  count_tmp +1;
end
end
sc_orb(sc_orb==1) = 0;
sc_vec = (sc_vec);
sc_orb = mod(sc_orb,1);
end
