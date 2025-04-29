function [WCC,WCCvec,HWan] = WannierCenter1D(WAVECAR_loop)
% WANNIERCENTER1D Calculate Wannier centers from 1D wavefunction data
%
% Syntax:
%   [WCC,WCCvec,HWan] = WannierCenter1D(WAVECAR_loop)
%
% Description:
%   Computes Wannier charge centers (WCCs) from wavefunction data along
%   a 1D k-loop using the Wilson loop method.
%
% Input Arguments:
%   WAVECAR_loop - Wavefunction data array (Norb x Nband x Nk)
%
% Output Arguments:
%   WCC - Wannier center values (eigenphases in [0,1))
%   WCCvec - Corresponding eigenvectors
%   HWan - Wilson loop matrix
%
% Example:
%   [wcc,vec] = WannierCenter1D(W);


%
% W(C)=M^{k0,k1}·...·M^{kn−1,kn,}
% M_{m, n}^{\mathbf{k}_{i}, \mathbf{k}_{j}}=\left\langle u_{m, \mathbf{k}_{i}} \mid u_{n, \mathbf{k}_{j}}\right\rangle
% M_{m, n}^{ki,kj} = <u_{m,ki}|u_{n,kj}>
% Take care!
% Bloch states : Psi_{n,k} = e^{i*k*r}*u_{n,k}
% So the u_{n,k} should be normalized phases instead of phi
HWan = eye(size(WAVECAR_loop,2));
Nloop = size(WAVECAR_loop,3);
for kj = 1:Nloop-1
    HWan = HWan* TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
end

%HWan = eye(size(WAVECAR_loop,2));
%Nloop = size(WAVECAR_loop,3);
%WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
%for kj = 1:Nloop-1
%    if kj == Nloop
%        F = TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,1));
%    else
%        F = TBkit.BerryConnection(WAVECAR_loop(:,:,kj),WAVECAR_loop(:,:,kj+1));
%    end
%    [U,~,V] = svd(F);
%    HWan = HWan*(U*V');
%    %Wan = Wan*F;
%end
%HWan = WAVECAR_loop(:,:,1)';
%Nloop = size(WAVECAR_loop,3);
%Wan = eye(size(WAVECAR_loop,1));
%WAVECAR_loop(:,:,Nloop+1) =WAVECAR_loop(:,:,1);
%for kj = 2:Nloop-1
%    Wan = Wan*(WAVECAR_loop(:,:,kj)*WAVECAR_loop(:,:,kj)');
%end
%HWan = HWan*Wan*WAVECAR_loop(:,:,end);%HWan*Wan*WAVECAR_loop(:,:,1)
%Wan
[WCCvec,WCCU] = eig(HWan);
%Ei = mod((angle(diag(WCCU))/(2*pi)),1);
Ei = mod(real(log(diag(WCCU))/(2*pi*1i)),1);
[WCCvec,WCC] = park.sorteig(Ei,WCCvec);
end