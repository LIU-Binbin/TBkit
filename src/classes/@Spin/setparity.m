function SpinObj = setparity(SpinObj, paritymat)
if isequal(size(SpinObj), size(paritymat))
for i = 1:numel(SpinObj)
SpinObj(i).parity = paritymat(i);
end
elseif length(paritymat) == 1
for i = 1:numel(SpinObj)
SpinObj(i).parity = paritymat;
end
end
end
