classdef Spin < HollowKnight
properties(Access = private, Dependent)
s
ms
end
properties(Hidden)
J
Jz
end
properties
parity = 1;
orientation = [0 0 1];
end
properties(Dependent, Hidden)
hollow;
s2_bar;
EigenVectors;
EigenVector;
end
methods
 SpinObj = Spin(S, Sz, coe, options)
end
methods
 A = uminus(A)
 C = minus(A, B, options)
 C = plus(A, B, options)
 C = innertimes(A, B)
 C = eq(A, B, options)
end
methods
 B = contractrow(A, options)
end
methods
 SpinObj = setparity(SpinObj, paritymat)
end
methods
 hollow = get.hollow(SpinObj)
 s = get.s(SpinObj)
 ms = get.ms(SpinObj)
 s2_bar = get.s2_bar(SpinObj)
 EigenVectors = get.EigenVectors(SpinObj)
end
methods
U = rotation(SpinObj, rotm, rightorleft, options)
U = rotation2(SpinObj, rotm, rightorleft, options)
end
methods
 Ak = rotateinner(A, abc, RightorLeft, immproper, conjugate, antisymmetry)
 Trmat = Tr(SpinObj)
 Invmat = ParityMat(SpinObj)
 WigerDmat = WignerD(SpinObj, abc, rightorleft)
 Matelement = WignerD_single(SpinObj1, SpinObj2, abc, rightorleft)
end
methods
 SpinObj = CG(Spinobj1, Spinobj2, options)
 SpinObj = TimeRerversal(SpinObj)
 SpinObj = Inversion(SpinObj)
 SpinObj = SzOper(Spinobj)
 SpinObj = SxOper(Spinobj)
 SpinObj = SyOper(Spinobj)
 SpinObj = SplusOper(Spinobj)
 SpinObj = SminusOper(Spinobj)
end
methods
 SzM = Sz(SpinObj,options)
 SyM = Sy(SpinObj,options)
 SxM = Sx(SpinObj,options)
 LzM = Lz(SpinObj,options)
 LxM = Lx(SpinObj,options)
 LyM = Ly(SpinObj,options)
 JzM = JZ(SpinObj,options)
 JxM = JX(SpinObj,options)
 JyM = JY(SpinObj,options)
 Splus = Splus(SpinObj,options)
 Sminus = Sminus(SpinObj,options)
 Lplus = Lplus(SpinObj,options)
 Lminus = Lminus(SpinObj,options)
 Jplus = Jplus(SpinObj,options)
 Jminus = Jminus(SpinObj,options)
end
methods
 SpinObj = BasisGen(Spinobj,force)
 SpinObj = PerM(Spinobj,permutation)
end
methods
 disp(SpinObj)
 Output = string(SpinObj)
 Output = pretty(SpinObj,options)
end
methods(Static)
end
methods(Static)
end
end
