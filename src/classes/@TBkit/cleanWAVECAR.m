function WAVECAR = cleanWAVECAR(WAVECAR,EIGENCAR,V,Accuracy)
%CLEANWAVECAR Smooth wavefunctions in degenerate subspaces
%
%   Syntax:
%       WAVECAR = cleanWAVECAR(WAVECAR,EIGENCAR,V,Accuracy)
%
%   Description:
%       Ensures smoothness of wavefunctions within degenerate subspaces
%       using perturbation theory.
%
%   Inputs:
%       WAVECAR  - Wavefunction array [orbitals×bands×k-points]
%       EIGENCAR - Eigenvalue array [bands×k-points]
%       V        - Perturbation potential
%       Accuracy - Degeneracy tolerance
%
%   Output:
%       WAVECAR - Smoothed wavefunctions
arguments
    WAVECAR;
    EIGENCAR;
    V;
    Accuracy = 1e-12;
end
kn = size(EIGENCAR,2);
%ki = 1;
%EIGEN_ki = EIGENCAR(:,ki);
%DegenPairBase = TBkit.checkDegen(EIGEN_ki,Accuracy);
%for ki = 2:kn
%    EIGEN_ki = EIGENCAR(:,ki);
%    DegenPair = TBkit.checkDegen(EIGEN_ki,Accuracy);
%    if ~isequal(DegenPairBase,DegenPair)
%        error('!');
%    end
%end
%WAVECAR = TBkit.smoothDegen2(WAVECAR);
for ki = 1:kn
    EIGEN_ki = EIGENCAR(:,ki);
    DegenPair = TBkit.checkDegen(EIGEN_ki,Accuracy);
    WAVE_ki = WAVECAR(:,:,ki);
    for i =1:size(DegenPair,1)
        % $\hat{H}^{(0)} \psi^{(1)}+\hat{H}^{\prime} \psi^{(0)}=E^{(0)} \psi^{(1)}+E^{(1)} \psi^{(0)}$
        WAVE_ki(:,DegenPair(i,1):DegenPair(i,2),:) = TBkit.smoothDegen(WAVE_ki(:,DegenPair(i,1):DegenPair(i,2)),V);
    end
    WAVECAR(:,:,ki) = WAVE_ki;
end
end
function WAVEFuncL = smoothDegen2(WAVEFuncL)
% Wave1 Todo
WAVEFuncLabs =abs(WAVEFuncL(:,1,:));
[~,maxWAVEFuncLabs] = max(WAVEFuncLabs,[],1);
[~,SelectOrb] = max(sum(maxWAVEFuncLabs,2));
for i = 1:size(WAVEFuncL,2)
    for k = 1:size(WAVEFuncL,3)
        WAVEFuncLtmp = exp(-1i*angle(WAVEFuncL(SelectOrb,1,i)))*WAVEFuncL(:,:,i);
        WAVEFuncL(:,:,i) = WAVEFuncLtmp;
    end
end