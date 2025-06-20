function [SYMCAR,OperObj] = Character_gen(OperObj,TBkitobj,klist,options)
arguments
    OperObj Oper;
    TBkitobj TBkit;
    klist = [0,0,0;0.5,0,0;0,0.5,0;0.5,0.5,0;0,0,0.5;0.5,0,0.5;0,0.5,0.5;0.5,0.5,0.5;];
    options.BasisFunction = [];
    options.generate_group = false;
    options.classify = true;
    options.Accuracy = 1e-6;
    options.sum = false;
end
Accuracy = round(log10(options.Accuracy));
if isempty(options.BasisFunction)
    BasisFunction = BasisFunc(TBkitobj);
else
    BasisFunction = options.BasisFunction;
end
if options.generate_group
    OperObj = OperObj.generate_group();
end
if options.classify
end
nOper = length(OperObj);
for i = 1:nOper
    if isnan(OperObj(i).U)
        OperObj(i) = OperObj(i).attachRm(TBkitobj.Rm);
        try
            OperObj(i).U = BasisFunction.rotation('Oper',OperObj(i));
        catch
            error('fail to generate Oper matrix!');
        end
    end
end
kn = size(klist,1);
switch class(TBkitobj)
    case 'HR'
        Basis_num = TBkitobj.WAN_NUM;
        [EIGENCAR,WAVECAR] = TBkitobj.EIGENCAR_gen('klist',klist,'convention','I','printmode',false);
    case {'Htrig','HK'}
        Basis_num = TBkitobj.Basis_num;
        [EIGENCAR,WAVECAR] = TBkitobj.EIGENCAR_gen('klist',klist*TBkitobj.Gk,'printmode',false);
end
VCAR = zeros(kn,Basis_num,nOper);
FactorL = zeros(kn,nOper);
orbL = TBkitobj.orbL;
for i = 1:nOper
    VCAR(:,:,i) = exp(1i*2*pi*(klist*OperObj(i).Rf-klist)*(orbL.'));
    FactorL(:,i) = exp(-1i*2*pi*(klist*OperObj(i).Rf)*OperObj(i).tf.');
end
SYMCAR = zeros(Basis_num,nOper,kn);
for i =1:kn
    for j = 1:nOper
        SYMCAR(:,j,i) = FactorL(i,j)*Oper.EIGEN_Kone(WAVECAR(:,:,i),OperObj(j).U,diag(VCAR(i,:,j)));
    end
end
SYMCAR = roundn(real(SYMCAR),Accuracy) + 1i*roundn(imag(SYMCAR),Accuracy);
if options.sum
end
end
