function mustBeEqualSize(a,b)
if ~(isequal(size(a),size(b)))
eid = 'Size:notEqual';
msg = 'Inputs must have same size.';
throwAsCaller(MException(eid,msg))
end
end
