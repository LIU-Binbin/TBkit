function BasisFunction = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
switch class(A.BFuncL)
case 'QnumL'
otherwise
BFuncL = rotateinner(A.BFuncL,abc,RightorLeft,immproper,conjugate,antisymmetry);
SpinL = rotateinner(A.spin,abc,RightorLeft,immproper,conjugate,antisymmetry);
BForb = BasisFunc.rotation_orb(A.BForb,R,t,options);
BasisFunction = BasisFunc();
end
end
