function EIGENCAR_SYM = EIGEN(OperObj,WAVECAR,klist,orbL,options)
arguments
OperObj Oper;
WAVECAR;
klist;
orbL;
options.Accuracy = 1e-6;
options.Obj = 'HR';
end
switch options.Obj
case 'HR'
VCAR = exp(1i*2*pi*(klist*OperObj.R-klist)*orbL.');
FactorL = exp(-1i*2*pi*klist*OperObj.R*OperObj.tf.');
case {'Htrig','HK'}
end
Accuracy = round(log10(options.Accuracy));
kn = size(WAVECAR,3);
NBANDS = size(WAVECAR,2);
EIGENCAR_SYM = zeros(NBANDS,kn);
for i =1:kn
EIGENCAR_SYM(:,i) = FactorL(i)*Oper.EIGEN_Kone(WAVECAR(:,:,i),OperObj.U,diag(VCAR(i,:)));
end
EIGENCAR_SYM = roundn(real(EIGENCAR_SYM),Accuracy) + 1i*roundn(imag(EIGENCAR_SYM),Accuracy);
end
