function EIGENCAR_Kone = EIGEN_Kone2(WAVECAR_one,U)
NBANDS = size(WAVECAR_one,2);
EIGENCAR_Kone = zeros(NBANDS,1);
for i = 1:NBANDS
EIGENCAR_Kone(i) = Oper.EIGEN_one(WAVECAR_one(:,i),U);
end
end
