function [pc_orb,pc_orbL_full,pc_elementL,pc_quantumL,orb_id_L,pc_orb_id_L,pc_orb_selectL] = unfold_orb(H_hr,Ns,Accuracy,orb_id_L)
arguments
H_hr HR;
Ns double = eye(3);
Accuracy double = 1e-6;
orb_id_L = [];
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
WANNUM = H_hr.WAN_NUM;
pc_orb=zeros(WANNUM,3);
translation_vector=zeros(WANNUM,3);
count_tmp = 1;
cur_sc_vec = [0 0 0];
for iorb = 1:WANNUM
pc_orb_tmp = roundn((orb_init(iorb,:)*Ns-cur_sc_vec),nAccuracy);
[pc_orb_incell,translation_vector(count_tmp,:)] = TBkit.translation_orb(pc_orb_tmp);
pc_orb(count_tmp,:) = pc_orb_incell;
count_tmp =  count_tmp +1;
end
pc_orbL_full = pc_orb;
PTL = [pc_orb,translation_vector];
[uniquePTL,uniqueAll,TrackingLater] = unique(PTL,'rows','stable');
translation_vector_mini = translation_vector(uniqueAll,:);
pc_orb_mini = pc_orb(uniqueAll,:);
[~,uniqueOrb] = uniquetol(pc_orb_mini,10^(nAccuracy+1),'ByRows',true);
[pc_orb_selectL,~] = ismember(PTL,uniquePTL(uniqueOrb,:),'rows');
pc_orb = pc_orbL_full(pc_orb_selectL,:);
pc_elementL = H_hr.elementL(pc_orb_selectL,:);
pc_quantumL = H_hr.quantumL(pc_orb_selectL,:);
pc_orb_id_L = 1:size(pc_orb,1);
if isempty(orb_id_L)
pc_peqL = [pc_orb,pc_elementL,pc_quantumL];
sc_peqL = [pc_orbL_full,H_hr.elementL,H_hr.quantumL];
else
pc_peqL = [pc_orb,orb_id_L(pc_orb_selectL).'];
sc_peqL = [pc_orbL_full,orb_id_L.'];
end
pc_peqL(isnan(pc_peqL)) = 0;
sc_peqL(isnan(sc_peqL)) = 0;
if size(unique(pc_peqL,"rows"),1) == size(pc_peqL,1)
else
warning('check duplicate orbital in primitive cell, the unfolding process is unreliable')
end
[~,sc_orb_selectL] = ismembertol(sc_peqL,pc_peqL,10^(nAccuracy+1),'ByRows',true);
orb_id_L = pc_orb_id_L(sc_orb_selectL);
end
