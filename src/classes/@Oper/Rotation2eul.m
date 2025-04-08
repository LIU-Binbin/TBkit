function eul = Rotation2eul( R, Rm )
arguments
R  {mustBeReal,Oper.mustBeSize(R,[3,3])}
Rm = [1 0 0;0 1 0;0 0 1];
end
if isequal(Rm,[1 0 0;0 1 0;0 0 1])
else
R = inv(Rm.' * R / Rm.');
end
if isa(R,'double')
try
eul = rotm2eul(R,'ZYZ');
catch
eul = zeros(1, 3, size(R,3), 'like', R);
nextAxis = [2, 3, 1, 2];
seqSettings.ZYZ = [3, 1, 1, 1];
setting = seqSettings.('ZYZ');
firstAxis = setting(1);
repetition = setting(2);
parity = setting(3);
movingFrame = setting(4);
i = firstAxis;
j = nextAxis(i+parity);
k = nextAxis(i-parity+1);
if repetition
sy = sqrt(R(i,j,:).*R(i,j,:) + R(i,k,:).*R(i,k,:));
singular = sy < 10 * eps(class(R));
eul = [atan2(R(i,j,:), R(i,k,:)), atan2(sy, R(i,i,:)), atan2(R(j,i,:), -R(k,i,:))];
numSingular = sum(singular,3);
assert(numSingular <= length(singular));
if numSingular > 0
eul(:,:,singular) = [atan2(-R(j,k,singular), R(j,j,singular)), ...
atan2(sy(:,:,singular), R(i,i,singular)), zeros(1,1,numSingular,'like',R)];
end
else
sy = sqrt(R(i,i,:).*R(i,i,:) + R(j,i,:).*R(j,i,:));
singular = sy < 10 * eps(class(R));
eul = [atan2(R(k,j,:), R(k,k,:)), atan2(-R(k,i,:), sy), atan2(R(j,i,:), R(i,i,:))];
numSingular = sum(singular,3);
assert(numSingular <= length(singular));
if numSingular > 0
eul(:,:,singular) = [atan2(-R(j,k,singular), R(j,j,singular)), ...
atan2(-R(k,i,singular), sy(:,:,singular)), zeros(1,1,numSingular,'like',R)];
end
end
if parity
eul = -eul;
end
if movingFrame
eul(:,[1,3],:)=eul(:,[3,1],:);
end
end
elseif isa(R,'sym')
firstAxis = 3;parity =1;movingFrame=1;
nextAxis = [2, 3, 1, 2];
i = firstAxis;
j = nextAxis(i+parity);
k = nextAxis(i-parity+1);
sy = simplify(sqrt(R(i,j,:).*R(i,j,:) + R(i,k,:).*R(i,k,:)));
eul = [atan2(R(i,j,:), R(i,k,:)), atan2(sy, R(i,i,:)), atan2(R(j,i,:), -R(k,i,:))];
singular = logical(sy ==sym(0));
numSingular = sum(singular,3);
assert(numSingular <= length(singular));
if numSingular > 0
eul(:,:,singular) = [angle(-R(j,k,singular)*1i+R(j,j,singular)), ...
angle(1i*sy(:,:,singular)+R(i,i,singular)), zeros(1,1,numSingular,'like',R)];
end
if parity
eul = -eul;
end
if movingFrame
eul(:,[1,3],:)=eul(:,[3,1],:);
end
eul = simplify(eul,'IgnoreAnalyticConstraints',true);
end
end
