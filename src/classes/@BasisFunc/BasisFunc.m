classdef BasisFunc < HollowKnight
properties
BForb = [];
BFnum =[];
BFuncL = [];
spin = [];
end
properties(Hidden,Dependent)
hollow;
end
properties(Access = private,Dependent)
FuncNum ;
BFclassL ;
end
properties(Hidden)
parity;
Rm = [1 0 0;0 1 0;0 0 1];
end
methods
 BasisFunction = BasisFunc(BFuncLOrigin,spin,BFnum,coe,orb,options)
end
methods(Static)
 BasisFunction = BasisFunction(TBkitobj)
end
methods
 hollow = get.hollow(BasisFunction)
 FuncNum = get.FuncNum(BasisFunction)
 BFclassL = get.BFclassL(BasisFunction)
end
methods
 dispAll(BasisFunction)
 disporbL(BasisFunction)
end
methods
 C = eq(A,B,options)
 C = innertimes(A,B,options)
end
methods
 B = contractrow(A,options)
end
methods(Static)
 [BFuncL] = introduce_coe(BFuncL,coeL)
 [coeL,BFuncL] = extract_coe(BFuncL,options)
end
methods
 U = rotation(A,Rc,Rf,tf,optionsConvection,optionsOper,optionsRm,options)
 Am = rotate(A,Rc,Rf,tf,rightorleft,optionsOper,options)
 A_Lj = rotaterow(A,Rc,Rf,tf,rightorleft,options)
 BasisFunction = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
end
methods(Static)
 [BFuncL,coeL] = rotation_func(BFuncL,R,t,options)
 spinL = rotation_spin(spin,R,t,options)
 BForb = rotation_orb(BForb,Rf,tf,options)
end
methods(Static)
end
end
