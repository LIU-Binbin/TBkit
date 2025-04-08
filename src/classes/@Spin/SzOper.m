function SpinObj = SzOper(Spinobj)
SpinObj = Spinobj;
for i = 1:numel(Spinobj)
SpinObj(i).coe = SpinObj(i).coe * SpinObj(i).coe * Spinobj(i).Jz;
end
SpinObj = SpinObj.contract();
end
