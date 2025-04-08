function SpinObj = SyOper(Spinobj)
SpinObj = 1/2i * (SplusOper(Spinobj) - SminusOper(Spinobj));
SpinObj = SpinObj.contract();
end
