function SpinObj = TimeRerversal(SpinObj)
for i = 1:numel(SpinObj)
SpinObj(i).coe = (-1)^(SpinObj(i).J - SpinObj(i).Jz) * SpinObj(i).coe;
SpinObj(i).Jz = -SpinObj(i).Jz;
end
end
