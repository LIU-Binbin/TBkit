function SpinObj = SplusOper(Spinobj)
SpinObj = Spinobj;
for i = 1:numel(Spinobj)
SpinObj(i).coe = SpinObj(i).coe * sqrt(Spinobj(i).s2_bar - Spinobj(i).Jz * (Spinobj(i).Jz + 1));
SpinObj(i).Jz = Spinobj(i).Jz + 1;
end
SpinObj = SpinObj.contract();
end
