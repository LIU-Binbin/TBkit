function SpinObj = BasisGen(Spinobj,force)
if nargin < 2
force = true;
end
if force
count = 0;
for Sz = Spinobj(1).J:-1:-Spinobj(1).J
count = count + 1;
SpinObj(count) = Spin(Spinobj(1).J, Sz, 'orientation', Spinobj(1).orientation);
end
end
end
