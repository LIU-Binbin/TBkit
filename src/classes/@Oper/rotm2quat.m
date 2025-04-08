function quat = rotm2quat( R )
if isa(R,'sym')
quat = sym(zeros(size(R,3), 4, 'like', R));
else
quat = zeros(size(R,3), 4, 'like', R);
end
K11 = R(1,1,:) - R(2,2,:) - R(3,3,:);
K12 = R(1,2,:) + R(2,1,:);
K13 = R(1,3,:) + R(3,1,:);
K14 = R(3,2,:) - R(2,3,:);
K22 = R(2,2,:) - R(1,1,:) - R(3,3,:);
K23 = R(2,3,:) + R(3,2,:);
K24 = R(1,3,:) - R(3,1,:);
K33 = R(3,3,:) - R(1,1,:) - R(2,2,:);
K34 = R(2,1,:) - R(1,2,:);
K44 = R(1,1,:) + R(2,2,:) + R(3,3,:);
K = [...
K11,    K12,    K13,    K14;
K12,    K22,    K23,    K24;
K13,    K23,    K33,    K34;
K14,    K24,    K34,    K44];
K = K ./ 3;
for i = 1:size(R,3)
[eigVec,eigVal] = eig(K(:,:,i));
eigVal = diag(eigVal);
[~,maxIdx] = max(real(eigVal));
quat(i,:) = real([eigVec(4,maxIdx) eigVec(1,maxIdx) eigVec(2,maxIdx) eigVec(3,maxIdx)]);
if quat(i,1) < 0
quat(i,:) = -quat(i,:);
end
end
end
