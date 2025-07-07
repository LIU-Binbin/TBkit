function dH_dk_xyz = dH_dk(H_hr, kpoint, options)
arguments
    H_hr HR;
    kpoint (1,3) double
    options.convention {mustBeMember(options.convention,{'I','II'})}= 'II';
end

vectorList = double(H_hr.vectorL);
vectorList_R = vectorList * H_hr.Rm;

Nbands = H_hr.Nbands;
dH_dk_R = zeros(Nbands,Nbands,3);

eikR = exp(1i*2*pi*vectorList*kpoint.'); % dim = R

for i = 1:3
    dH_dk_R(:,:,i) = 1i * tensorprod(H_hr.HnumL, eikR.*vectorList_R(:,i), 3, 1);
end
%% not tested!
if options.convention == 'I'   
    tji_mat_cart = H_hr.tjmti{1};
    tji_mat_frac = H_hr.tjmti{2};
    eikr = exp(1i*2*pi*(tji_mat_frac*kpoint.')); % dim=ij

    % (H_R*exp(ik(R+r)))' = (H_R*exp(ikR)*ikR) * exp(ikr) + (H_R*exp(ikR)) * ikr*exp(ikr)
    Hk = tensorprod(H_hr.HnumL, eikR, 3, 1); % dim=ij
    dH_dk_xyz = dH_dk_R .* eikr + Hk .* (1i*tji_mat_cart) .* eikr;
else
    dH_dk_xyz = dH_dk_R;
end
end