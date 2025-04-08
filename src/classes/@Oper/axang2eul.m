function eul = axang2eul(axang)
R = Oper.axang2rotm( axang );
eul = Oper.Rotation2eul(R);
end
