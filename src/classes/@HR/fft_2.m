function [W,D,dH_dk_R,dH_dk_dk_R] = fft_2(H_hr, klist_cart, precomputed)
arguments
    H_hr HR
    klist_cart % cart
    precomputed = [] ; % 可选缓存
end

% try
%     [W, D, dH_dk_R] = fft_mex(H_hr.HnumL, double(H_hr.vectorL(:,1:H_hr.Dim))* H_hr.Rm, klist);
% catch
%     HnumList = H_hr.HnumL ;
%     vectorList = double(H_hr.vectorL(:,1:H_hr.Dim)) ;
%     vectorList_R = vectorList * H_hr.Rm;
%     FactorList = exp(1i*vectorList_R*klist.');
%     Hout = tensorprod(HnumList, FactorList, 3, 1);
%     [W,D]= eig((Hout+Hout')/2,'vector');
%     for i = 1:3
%         dH_dk_R(:,:,i) = 1i * tensorprod(HnumList , FactorList.*vectorList_R(:,i), 3, 1); % only for one kpoint?
%     end
% end
% % return;
persistent HnumList vectorList_R 
if isempty(HnumList)
    HnumList = H_hr.HnumL;
    vectorList = double(H_hr.vectorL(:,1:H_hr.Dim));  % M×3
    vectorList_R = vectorList * H_hr.Rm;  % M×3 real-space R
end


% 计算因子：exp(i * R · k)
phase = exp(1i * (vectorList_R * klist_cart(:)));  % M×1
Hk = tensorprod(HnumList, phase, 3, 1);
Hk = (Hk + Hk') / 2;  % Hermitian 修正
[W,D] = eig(Hk, 'vector');
for i = 1:3
    dH_dk_R(:,:,i) = 1i * tensorprod(HnumList , phase.*vectorList_R(:,i), 3, 1);
    for j = 1:3
        dH_dk_dk_R(:,:,i,j) = 1i * tensorprod(HnumList , phase.*vectorList_R(:,i).*vectorList_R(:,j), 3, 1);
    end
end


end
