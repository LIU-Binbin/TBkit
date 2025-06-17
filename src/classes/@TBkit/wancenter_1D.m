function [BF,BF_WAVE,Wan] = wancenter_1D(WAVECAR_loop,mode)
% WANCENTER_1D Calculate Wannier centers for 1D systems
%
% Syntax:
%   [BF,BF_WAVE,Wan] = wancenter_1D(WAVECAR_loop)
%   [BF,BF_WAVE,Wan] = wancenter_1D(WAVECAR_loop,mode)
%
% Description:
%   Computes Wannier centers (Wannier charge centers) from wavefunction data
%   along a 1D k-path using different formulations of the Wilson loop.
%
% Input Arguments:
%   WAVECAR_loop - Wavefunction data array (Norb x Nband x Nk)
%   mode - Calculation mode (optional, 'nested' for alternative formulation)
%
% Output Arguments:
%   BF - Wannier center values (eigenphases)
%   BF_WAVE - Corresponding eigenvectors
%   Wan - Wilson loop matrix
%
% Example:
%   [wcc,vec] = wancenter_1D(W);
formula = 2;
if nargin < 2
    switch formula
        case 1
            Nband = size(WAVECAR_loop,2);
            Wan = eye(Nband);
            Nloop = size(WAVECAR_loop,3);
            %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
            for k =1:Nloop-1
                WanTmp = zeros(Nband) ;
                for i = 1:Nband
                    for j = 1:Nband
                        WanTmp(i,j) = WAVECAR_loop(:,j,k+1).'*conj(WAVECAR_loop(:,i,k));
                    end
                end
                Wan = Wan*WanTmp;
            end
            %pause(1);
        case 2
            Wan = eye(size(WAVECAR_loop,2));
            Nloop = size(WAVECAR_loop,3);
            %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
            for kj = 1:Nloop-1
                Wan = Wan* TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
            end
            %Wan = Wan;%* TBkit.BerryConnection(WAVECAR_loop(:,:,end),WAVECAR_loop(:,:,1));
        case 3
            Wan = eye(size(WAVECAR_loop,2));
            Nloop = size(WAVECAR_loop,3);
            %WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
            for kj = 1:Nloop-1
                if kj == Nloop
                    F = TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,1));
                else
                    F = TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
                end
                [U,~,V] = svd(F);
                Wan = Wan*(U*V');
                %Wan = Wan*F;
            end
    end
    switch formula
        case 1
            [BF_WAVE,Ei] = eig(Wan);
            Ei = angle((Ei));
            [BF_WAVE,Ei] = sorteig(Ei,BF_WAVE);
            BF=diag(Ei);
        case {2,3}
            [BF_WAVE,Ei] = eig(Wan);
            Ei = angle((Ei));
            [BF_WAVE,Ei] = sorteig(Ei,BF_WAVE);
            BF=diag(Ei);
    end

elseif strcmp(mode,'nested')
    %                 if size(site_weight,2) ==size(WAVECAR_loop,3)
    %                     Wan = (WAVECAR_loop(:,:,1).*site_weight(:,1))';
    %                     for kj = 2:size(WAVECAR_loop,3)
    %                         Wan = Wan* (WAVECAR_loop(:,:,kj).*site_weight(:,kj))*(WAVECAR_loop(:,:,kj).*site_weight(:,kj))';
    %                     end
    %                     Wan = Wan* WAVECAR_loop(:,:,1).*site_weight(:,1);
    %                 else
    %                     Wan = (WAVECAR_loop(:,:,1).*site_weight(:,1))';
    %                     for kj = 2:size(WAVECAR_loop,3)
    %                         Wan = Wan* (WAVECAR_loop(:,:,kj).*site_weight(:,1))*(WAVECAR_loop(:,:,kj).*site_weight(:,1))';
    %                     end
    %                     Wan = Wan* (WAVECAR_loop(:,:,1).*site_weight(:,1));
    %                 end
    Wan = WAVECAR_loop(:,:,1)';
    for kj = 2:size(WAVECAR_loop,3)
        Wan = Wan* WAVECAR_loop(:,:,kj)*WAVECAR_loop(:,:,kj)';
    end
    Wan = Wan* WAVECAR_loop(:,:,1);
    [BF_WAVE,Ei] = eig(Wan);
    BF=angle(sum(diag(Ei)));
end