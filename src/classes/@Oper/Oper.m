classdef Oper < group
properties
R =[];
t = [0 0 0];
conjugate = false;
antisymmetry = false;
Rf = [];
tf = [0 0 0];
end
properties (GetAccess = protected,Hidden = true)
continuous = false;
strict_eq = false;
end
properties (GetAccess = protected)
end
methods
 SymOper = Oper(R,U,t,options)
end
methods(Static)
 SymOper = identity(dim, shape,propArgs)
 SymOper = time_reversal(realspace_dim, U, spin,propArgs)
 SymOper = particle_hole(realspace_dim, U,propArgs)
 SymOper = chiral(realspace_dim, U,propArgs)
 SymOper = inversion(realspace_dim, U,quantumL,propArgs)
 SymOper = rotation(angle, axis, inversion, U, spin,options,propArgs)
 SymOper = spaceRotation(angle, axis,t, inversion, U, spin,options,propArgs)
 SymOper = C3z(realspace_dim, inversion, U, spin)
 SymOper = C4z(realspace_dim, inversion, U, spin)
 SymOper = C6z(realspace_dim, inversion, U, spin)
 SymOper = mirror(axis, U, spin,options,propArgs)
 Mx
 My
 Mz
 group = square(tr, ph, generators, spin,options,propArgs)
 group = cubic(tr, ph, generators, spin,options,propArgs)
 group = hexagonal_2D(tr, ph, generators, spin,options, propArgs)
 group = hexagonal(tr, ph, generators, spin,options,propArgs)
end
methods
 [TrueOrNot,result] = commute(SymOper1,SymOper2)
end
methods
 SymOper = Ugen(SymOper,Basis,options)
end
methods
 [SYMCAR,OperObj] = Character_gen(OperObj,TBkitobj,klist,options)
 EIGENCAR_SYM = EIGEN(OperObj,WAVECAR,klist,orbL,options)
end
methods(Static)
 EIGENCAR_Kone = EIGEN_Kone(WAVECAR_one,D,V)
 EIGENCAR_Kone = EIGEN_Kone2(WAVECAR_one,U)
 EIGEN = EIGEN_one(WAVEFUNC,U)
end
methods(Static)
 U = Umat(SymOper1,SymOper2,options)
 P = similarity(SymOper1,SymOper2,options)
end
methods
 str = disp(SymOper,options)
 [SymMat,SymR] = sym(SymOper)
 SymOper = E(SymOper)
end
methods
 [basic_eq,U_eq]=eq(SymOper1,SymOper2)
 result = lt(SymOper1,SymOper2)
 SymOper_out = times(SymOper1,SymOper2)
 SymOper = inv(SymOper)
end
methods
 OperObj = attachRm(OperObj,Rm)
end
methods
 SymOper_Str = string(SymOper)
 SymOper_latex= latex(SymOper)
 SymOper_pretty= repr_pretty(SymOper,Cycle)
 SymOper_latex= repr_latex(SymOper)
 name= pretty(SymOper,options)
end
methods
 symmetry_from_permutation()
end
methods (Access= protected)
end
methods (Static)
 S_mat = spin_matrices(s, include_0)
 J_mat = spin_matrices_from_orb(quantumL,include_0,strict)
 L_mat = L_matrices(d, l)
 U = spin_rotation(n, s, roundint)
 phi = braket(phi1,Oper,phi2)
end
methods(Static)
 angle = name_angle(theta, Latex)
 str = mat2latex(mat,accuracy)
 str = num2latex(num,accuracy)
 C = isclose(A,B)
 C = allclose(A,B)
 n = round_axis(n)
 M = tensordot_naive(A,B,sizeLast)
 RotationMat = nThetad2Rotation(thetad,n)
 RotationMat = nTheta2Rotation(theta,n)
 [n,theta]= Rotation2nTheta(R,Rm)
 [n,thetad]= Rotation2nThetad(R,Rm)
 Rf = Rc2Rf(Rc,Rm)
 Rft = Rct2Rft(Rct,Rm)
 Rc = Rf2Rc(Rf,Rm)
 Rct = Rft2Rct(Rft,Rm)
 [prop,Coeff] = prop_to_id(A)
 mustBeOfClass(input,className)
 mustBeEqualSize(a,b)
 mustBeSize(a,b)
 mustBeDims(input,numDims)
 mustBeHalfInteger(A)
 [Asort,Usort] = sorteig(U,A)
end
methods(Static)
 eul = Rotation2eul( R, Rm )
 R = axang2rotm( axang )
 eul = axang2eul(axang)
 quat = rotm2quat( R )
end
end
