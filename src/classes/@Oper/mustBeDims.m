function mustBeDims(input,numDims)
if ~isequal(length(size(input)),numDims)
eid = 'Size:wrongDimensions';
msg = ['Input must have dimensions: ',num2str(numDims)];
throwAsCaller(MException(eid,msg))
end
end
