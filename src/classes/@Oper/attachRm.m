function OperObj = attachRm(OperObj,Rm)
for i =1:numel(OperObj)
OperObj(i).Rf = double(Oper.Rc2Rf(OperObj(i).R,Rm));
OperObj(i).tf = OperObj(i).t;
end
end
