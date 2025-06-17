function WigerDmat = WignerD(SpinObj, abc, rightorleft)
arguments
    SpinObj Spin;
    abc;
    rightorleft = 'right';
end
nSpin = size(SpinObj, 1);
WigerDmat = sym(zeros(nSpin));
for i = 1:nSpin
    for j = 1:nSpin
        Matelement = WignerD_single(SpinObj(i,:), SpinObj(j,:), abc, rightorleft);
        if Matelement ~= sym(0)
            WigerDmat(i, j) = Matelement;
        end
    end
end
end
