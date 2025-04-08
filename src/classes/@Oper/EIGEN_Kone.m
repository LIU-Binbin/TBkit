function EIGENCAR_Kone = EIGEN_Kone(WAVECAR_one,D,V)
EIGENCAR_Kone = diag(WAVECAR_one'*V*D*WAVECAR_one);
end
