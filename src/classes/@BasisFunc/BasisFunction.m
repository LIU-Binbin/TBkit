function BasisFunction = BasisFunction(TBkitobj)
switch class(TBkitobj)
case {'TBkit','HR','Htrig','HK'}
[BFuncLOrigin,S,SzL] = Qnum.QnumL(TBkitobj);
spinL = Spin(S,SzL);
BasisFunction = BasisFunc(BFuncLOrigin,spinL,1,1,TBkitobj.orbL);
for i = 1:numel(BasisFunction)
BasisFunction(i).Rm = TBkitobj.Rm;
end
end
end
