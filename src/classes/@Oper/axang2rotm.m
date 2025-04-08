function R = axang2rotm( axang )
if isa(axang,'double')
v = normalize(axang(:,1:3),'scale');
else
v =(axang(:,1:3));
end
numInputs = size(axang,1);
theta = zeros(1,1,numInputs,class(axang));
theta(1,1,:) = axang(:,4);
cth = cos(theta);
sth = sin(theta);
vth = (1 - cth);
vx = zeros(1,1,numInputs,'like',axang);
vy = vx;
vz = vx;
vx(1,1,:) = v(:,1);
vy(1,1,:) = v(:,2);
vz(1,1,:) = v(:,3);
tempR = cat(1, vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth, ...
vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth, ...
vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth);
R = reshape(tempR, [3, 3, length(vx)]);
R = permute(R, [2 1 3]);
if isa(R, 'sym')
R = simplify(R);
end
end
