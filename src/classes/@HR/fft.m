function [W,D,dH_dk_R] = fft(H_hr, klist_cart, precomputed)
arguments
    H_hr HR
    klist_cart % cart
    precomputed = [] ; % 可选缓存
end
persistent HnumList vectorList_R NBANDS NRPTS
if isempty(HnumList)
    HnumList = H_hr.HnumL;
    vectorList = double(H_hr.vectorL(:,1:H_hr.Dim));  % M×3
    vectorList_R = vectorList * H_hr.Rm;  % M×3 real-space R
    NBANDS = H_hr.WAN_NUM;
    NRPTS = H_hr.NRPTS;
end


% 计算因子：exp(i * R · k)
phase = exp(1i * (vectorList_R * klist_cart(:)));  % M×1


% Hk = zeros(NBANDS, NBANDS, 'like', phase);  % 初始化哈密顿量
% for i = 1:NRPTS
%     Hk = Hk + HnumList(:,:,i) * phase(i);
% end

Hk = tensorprod(HnumList, phase, 3, 1);

Hk = (Hk + Hk') / 2;  % Hermitian 修正
[W,D] = eig(Hk, 'vector');


for i = 1:3
    dH_dk_R(:,:,i) = 1i * tensorprod(HnumList , phase.*vectorList_R(:,i), 3, 1);
    %dH_dk_R(:,:,i) = dH;
    %dH_dk_R(:,:,i) = (dH + dH') / 2;  % Hermitian修正 ?
end

% dH_dk_R = zeros(NBANDS, NBANDS, 3);
% for dir = 1:3  % x,y,z方向
%     dH = zeros(NBANDS, NBANDS);
%     for i = 1:NRPTS
%         % 导数项: i * R_dir * H(R) * exp(iR·k)
%         dH = dH + HnumList(:,:,i) * (1i * vectorList_R(i, dir) * phase(i));
%     end
%     %dH_dk_R(:,:,dir) = dH;
%     dH_dk_R(:,:,dir) = (dH + dH') / 2;  % Hermitian修正
% end

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