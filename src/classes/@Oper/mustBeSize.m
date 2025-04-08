function mustBeSize(a,b)
if sum(~(size(a)==b))
eid = 'Size:notRequired';
msg = 'Inputs must have size of '+mat2str(b);
throwAsCaller(MException(eid,msg))
end
end
