function Trmat = Tr(SpinObj)
Trmat = InnerProduct(SpinObj, TimeRerversal(SpinObj));
end
