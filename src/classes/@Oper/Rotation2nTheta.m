function [n,theta]= Rotation2nTheta(R,Rm)
arguments
R  {Oper.mustBeSize(R,[3,3])}
Rm = [1 0 0;0 1 0;0 0 1];
end
if isequal(Rm,[1 0 0;0 1 0;0 0 1])
else
R = inv(Rm.' * R / Rm.');
end
if Oper.isclose(det(R),-1)
R = -R;
disp('det = -1');
else
R = R/det(R);
end
if Oper.allclose(R,eye(3))
theta = 0;
n = [1, 0 ,0];
return
end
L_mat = Oper.L_matrices();
n = real(-1i*[trace(L_mat(:,:,1)*R),trace(L_mat(:,:,2)*R),trace(L_mat(:,:,3)*R)]);
absn = norm(n) * sign(sum(n));
if Oper.isclose(absn,0)
if isa(R,'sym')
[val, vec] = eig(simplify(R));
else
[val, vec] = eig((R));
end
[val,vec] =  Oper.sorteig(vec,val);
vec = diag(vec);
if Oper.isclose(vec(end),1)
n = val(:, end).';
n = n/ sign(sum(n));
end
else
n = n/absn;
end
try
theta = real(acos(complex((1/2)*(R(1,1,:)+R(2,2,:)+R(3,3,:)-1))));
catch
theta = real(acos(((1/2)*(R(1,1,:)+R(2,2,:)+R(3,3,:)-1))));
end
if isa(n, 'sym')
n = simplify(n,'IgnoreAnalyticConstraints',true);
theta = simplify(theta,'IgnoreAnalyticConstraints',true);
end
end
