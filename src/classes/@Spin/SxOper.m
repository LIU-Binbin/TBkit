function SpinObj = SxOper(Spinobj)
SpinObj = 1/2 * (SplusOper(Spinobj) + SminusOper(Spinobj));
SpinObj = SpinObj.contract();
end
