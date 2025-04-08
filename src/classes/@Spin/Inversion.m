function SpinObj = Inversion(SpinObj)
for i = 1:numel(SpinObj)
SpinObj(i).coe = SpinObj(i).coe * SpinObj(i).parity;
end
end
