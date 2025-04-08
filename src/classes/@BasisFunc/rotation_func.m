function [BFuncL,coeL] = rotation_func(BFuncL,R,t,options)
arguments
BFuncL
R
t
options.conjugate logical = false ;
options.Rm  = [1 0 0;0 1 0;0 0 1] ;
end
optionsCell = namedargs2cell(options);
for i =numel(BFuncL)
if isa(BFuncL{i},'sym')
else
BFuncL{i} = rotate(BFuncL{i},R,t,optionsCell{:});
end
end
end
