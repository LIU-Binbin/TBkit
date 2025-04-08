function Invmat = ParityMat(SpinObj)
Invmat = InnerProduct(SpinObj, Inversion(SpinObj));
end
