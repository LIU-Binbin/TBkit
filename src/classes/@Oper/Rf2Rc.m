function Rc = Rf2Rc(Rf,Rm)
rotation_cart_inv = Rm.' * rotation / Rm.';
Rc = inv(rotation_cart_inv );
end
