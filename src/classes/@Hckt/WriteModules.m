function WriteModules(Modules,fid,magnitude,Lprefix)
arguments
Modules char;
fid;
magnitude = 'p';
Lprefix = '1';
end
switch magnitude
case 'p'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = '';
case 'n'
Cmagnitude = 'n';
Lmagnitude = 'u';
Rmagnitude = '';
case 'u'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = '';
case 'm'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = '';
end
fprintf(fid,"* ---------------\n");
switch Modules
case 'Basis'
fprintf(fid,"* BasisC3_origin \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND VarL0=%f%s InitV=0V R_L=1u \n",Lprefix,Lmagnitude);
fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
fprintf(fid,"L3 n3 n2 VarL0 R=R_L \n");
fprintf(fid,".ends BasisC3_origin\n");
case 'Basis_SOI'
fprintf(fid,"* BasisC3_origin \n");
fprintf(fid,"*\n");
fprintf(fid,"* n1 n2 n3 n1_prime n2_prime n3_prime nR1_prime nR2_prime nR3_prime  \n");
fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND " + ...
"VarL0=1.8
fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
fprintf(fid,"L3 n3 n1 VarL0 R=R_L \n");
fprintf(fid,"X1_1prime n1 n1_prime TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime n2 n2_prime TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime n3 n3_prime TOGND VoltageFollower \n");
fprintf(fid,"X1prime n2_prime n1_prime n3_prime NR1_prime TOGND AdderSubtractor \n");
fprintf(fid,"X2prime n3_prime n2_prime n1_prime NR2_prime TOGND AdderSubtractor \n");
fprintf(fid,"X3prime n1_prime n3_prime n2_prime NR3_prime TOGND AdderSubtractor \n");
fprintf(fid,"R1prime n2 NR1_prime VarRh \n");
fprintf(fid,"R2prime n3 NR2_prime VarRh \n");
fprintf(fid,"R3prime n1 NR3_prime VarRh \n");
fprintf(fid,"C1 n1 TOGND  VarCg \n");
fprintf(fid,"C2 n2 TOGND  VarCg \n");
fprintf(fid,"C3 n3 TOGND  VarCg \n");
fprintf(fid,".ends BasisC3_origin\n");
case 'Basis_Chern1'
fprintf(fid,"* Basis_M1_module \n");
fprintf(fid,"*\n");
fprintf(fid,"* n1 n2 n3 n1_prime n2_prime n3_prime \n");
fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND " + ...
"VarRm=30
fprintf(fid,"R1 n2 n1_prime VarRm\n");
fprintf(fid,"R2 n3 n2_prime VarRm\n");
fprintf(fid,"R3 n1 n3_prime VarRm\n");
fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
fprintf(fid,"L3 n3 n1 VarL0 R=R_L \n");
fprintf(fid,"X1_prime n1 n1_prime TOGND InvertingOpAmp \n");
fprintf(fid,"X2_prime n2 n2_prime TOGND InvertingOpAmp \n");
fprintf(fid,"X3_prime n3 n3_prime TOGND InvertingOpAmp \n");
fprintf(fid,".ends BasisC3_origin\n");
case 'Basis_Chern2'
fprintf(fid,"* Basis_M2_module \n");
fprintf(fid,"*\n");
fprintf(fid,"* n1 n2 n3 \n");
fprintf(fid,".SubCkt BasisC3_origin n1 n2 n3 TOGND " + ...
"Var2Rm=30
fprintf(fid,"R1 n1 TOGND Var2Rm\n");
fprintf(fid,"R2 n2 TOGND Var2Rm\n");
fprintf(fid,"R3 n3 TOGND Var2Rm\n");
fprintf(fid,"L1 n1 n2 VarL0 R=R_L \n");
fprintf(fid,"L2 n2 n3 VarL0 R=R_L \n");
fprintf(fid,"L3 n3 n1 VarL0 R=R_L \n");
fprintf(fid,"C1 n1 n1_prime Var3Cm\n");
fprintf(fid,"C2 n2 n2_prime Var3Cm\n");
fprintf(fid,"C3 n3 n3_prime Var3Cm\n");
fprintf(fid,"X1_prime n2 n3_prime TOGND IntegratorOpAmp \n");
fprintf(fid,"X2_prime n3 n1_prime TOGND IntegratorOpAmp \n");
fprintf(fid,"X3_prime n1 n2_prime TOGND IntegratorOpAmp \n");
fprintf(fid,".ends BasisC3_origin\n");
case '+sigma_0'
fprintf(fid,"* +sigma_0 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1 R_n1  VarC0 IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n2 VarC0 IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n3 VarC0 IC=InitV\n");
fprintf(fid,".ends PlusSigma0\n");
case '-sigma_0'
fprintf(fid,"* -sigma_0 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1 R_n2 VarC0 IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n1 VarC0 IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n1 VarC0 IC=InitV\n");
fprintf(fid,"C4 L_n1 R_n3 VarC0 IC=InitV\n");
fprintf(fid,"C5 L_n2 R_n3 VarC0 IC=InitV\n");
fprintf(fid,"C6 L_n3 R_n2 VarC0 IC=InitV\n");
fprintf(fid,".ends MinusSigma0\n");
case '+isigma_0'
fprintf(fid,"* +isigma_0 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1  R_n2 VarR0\n");
fprintf(fid,"R2 L_n2 R_n1  VarR0\n");
fprintf(fid,"R3 L_n3 R_n1  VarR0\n");
fprintf(fid,"R4 L_n1  R_n3 VarR0\n");
fprintf(fid,"R5 L_n2 R_n3 VarR0\n");
fprintf(fid,"R6 L_n3 R_n2 VarR0\n");
fprintf(fid,".ends PlusiSigma0\n");
case '-isigma_0'
fprintf(fid,"* -isigma_0 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiSigma0 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1  R_n1  VarR0\n");
fprintf(fid,"R2 L_n2 R_n2 VarR0\n");
fprintf(fid,"R3 L_n3 R_n3 VarR0\n");
fprintf(fid,".ends MinusiSigma0\n");
case '+sigma_1'
fprintf(fid,"* +sigma_1 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=6.8
fprintf(fid,"C1 L_n1 R_n2 VarC0 IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n1 VarC0 IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n3 VarC0 IC=InitV\n");
fprintf(fid,".ends PlusSigma1\n");
case '-sigma_1'
fprintf(fid,"* -sigma_1 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=6.8
fprintf(fid,"C1 L_n1 R_n1 VarC0 IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n2 VarC0 IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n1 VarC0 IC=InitV\n");
fprintf(fid,"C4 L_n1 R_n3 VarC0 IC=InitV\n");
fprintf(fid,"C5 L_n2 R_n3 VarC0 IC=InitV\n");
fprintf(fid,"C6 L_n3 R_n2 VarC0 IC=InitV\n");
fprintf(fid,".ends MinusSigam1\n");
case '+isigma_1'
fprintf(fid,"* +isigma_1 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1  R_n1  VarR0\n");
fprintf(fid,"R2 L_n2 R_n2 VarR0\n");
fprintf(fid,"R3 L_n3 R_n1  VarR0\n");
fprintf(fid,"R4 L_n1  R_n3 VarR0\n");
fprintf(fid,"R5 L_n2 R_n3 VarR0\n");
fprintf(fid,"R6 L_n3 R_n2 VarR0\n");
fprintf(fid,".ends PlusiSigma1\n");
case '-isigma_1'
fprintf(fid,"* -isigma_1 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiSigma1 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1  R_n2 VarR0\n");
fprintf(fid,"R2 L_n2 R_n1  VarR0\n");
fprintf(fid,"R3 L_n3 R_n3 VarR0\n");
fprintf(fid,".ends MinusiSigma1\n");
case '+gen3sigma_2'
fprintf(fid,"* +gen3sigma_2 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1  R_n2 VarC0  IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n1  VarC0  IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n3 VarC0  IC=InitV\n");
fprintf(fid,"C4 L_n1  R_n3 Var2C0 IC=InitV\n");
fprintf(fid,"C5 L_n2 R_n2 Var2C0 IC=InitV\n");
fprintf(fid,"C6 L_n3 R_n1  Var2C0 IC=InitV\n");
fprintf(fid,".ends PlusGen3Sigma2\n");
case '-gen3sigma_2'
fprintf(fid,"* -gen3sigma_2 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1  R_n2 VarC0  IC=InitV\n");
fprintf(fid,"C2 L_n2 R_n1  VarC0  IC=InitV\n");
fprintf(fid,"C3 L_n3 R_n3 VarC0  IC=InitV\n");
fprintf(fid,"C4 L_n1  R_n1  Var2C0 IC=InitV\n");
fprintf(fid,"C5 L_n2 R_n3 Var2C0 IC=InitV\n");
fprintf(fid,"C6 L_n3 R_n2 Var2C0 IC=InitV\n");
fprintf(fid,".ends MinusGen3Sigma2\n");
case '+igen3sigma_2'
fprintf(fid,"* +igen3sigma_2 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1 R_n2  VarR0  \n");
fprintf(fid,"R2 L_n2 R_n1  VarR0  \n");
fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
fprintf(fid,"R4 L_n1 R_n1   VarR0_2\n");
fprintf(fid,"R5 L_n2 R_n3 VarR0_2\n");
fprintf(fid,"R6 L_n3 R_n2 VarR0_2\n");
fprintf(fid,".ends PlusiGen3Sigma2\n");
case '-igen3sigma_2'
fprintf(fid,"* -igen3sigma_2 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiGen3Sigma2 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1 R_n2  VarR0  \n");
fprintf(fid,"R2 L_n2 R_n1  VarR0  \n");
fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
fprintf(fid,"R4 L_n1 R_n3  VarR0_2\n");
fprintf(fid,"R5 L_n2 R_n2 VarR0_2\n");
fprintf(fid,"R6 L_n3 R_n1  VarR0_2\n");
fprintf(fid,".ends MinusiGen3Sigma2\n");
case '+gen3sigma_3'
fprintf(fid,"* +gen3sigma_3 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1  R_n1  VarR0  \n");
fprintf(fid,"R2 L_n2 R_n2 VarR0  \n");
fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
fprintf(fid,"R4 L_n1  R_n2 VarR0_2\n");
fprintf(fid,"R5 L_n2 R_n3 VarR0_2\n");
fprintf(fid,"R6 L_n3 R_n1  VarR0_2\n");
fprintf(fid,".ends PlusGen3Sigma3\n");
case '-gen3sigma_3'
fprintf(fid,"* -gen3sigma_3 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1 R_n1   VarR0  \n");
fprintf(fid,"R2 L_n2 R_n2 VarR0  \n");
fprintf(fid,"R3 L_n3 R_n3 VarR0  \n");
fprintf(fid,"R4 L_n1 R_n3  VarR0_2\n");
fprintf(fid,"R5 L_n2 R_n1  VarR0_2\n");
fprintf(fid,"R6 L_n3 R_n2 VarR0_2\n");
fprintf(fid,".ends MinusGen3Sigma3\n");
case '+igen3sigma_3'
fprintf(fid,"* +igen3sigma_3 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
fprintf(fid,"C4 L_n1 R_n2 Var2C0\n");
fprintf(fid,"C5 L_n2 R_n3 Var2C0\n");
fprintf(fid,"C6 L_n3 R_n1 Var2C0\n");
fprintf(fid,".ends PlusiGen3Sigma3\n");
case '-igen3sigma_3'
fprintf(fid,"* -igen3sigma_3 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1 R_n2 R_n3 TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1 R_n1 VarC0 \n");
fprintf(fid,"C2 L_n2 R_n2 VarC0 \n");
fprintf(fid,"C3 L_n3 R_n3 VarC0 \n");
fprintf(fid,"C4 L_n1 R_n3 Var2C0\n");
fprintf(fid,"C5 L_n2 R_n1 Var2C0\n");
fprintf(fid,"C6 L_n3 R_n2 Var2C0\n");
fprintf(fid,".ends MinusiGen3Sigma3\n");
case '+isigma_1_SOI'
fprintf(fid,"* +isigma_1_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiSigma1 L_n1_prime L_n2_prime L_n3_prime R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1 R_n1_prime VarR0\n");
fprintf(fid,"R2 L_n2 R_n2_prime VarR0\n");
fprintf(fid,"R3 L_n3 R_n1_prime VarR0\n");
fprintf(fid,"R4 L_n1 R_n3_prime VarR0\n");
fprintf(fid,"R5 L_n2 R_n3_prime VarR0\n");
fprintf(fid,"R6 L_n3 R_n2_prime VarR0\n");
fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,"R1_ L_n1_prime R_n2 VarR0\n");
fprintf(fid,"R2_ L_n2_prime R_n1 VarR0\n");
fprintf(fid,"R3_ L_n3_prime R_n3 VarR0\n");
fprintf(fid,"X1_1prime_ R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends PlusiSigma1\n");
case '-isigma_1_SOI'
fprintf(fid,"* -isigma_1_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiSigma1 L_n1_prime L_n2_prime L_n3_prime R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1_prime R_n2 VarR0\n");
fprintf(fid,"R2 L_n2_prime R_n1 VarR0\n");
fprintf(fid,"R3 L_n3_prime R_n3 VarR0\n");
fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,"R1_ L_n1 R_n1_prime VarR0\n");
fprintf(fid,"R2_ L_n2 R_n2_prime VarR0\n");
fprintf(fid,"R3_ L_n3 R_n1_prime VarR0\n");
fprintf(fid,"R4_ L_n1 R_n3_prime VarR0\n");
fprintf(fid,"R5_ L_n2 R_n3_prime VarR0\n");
fprintf(fid,"R6_ L_n3 R_n2_prime VarR0\n");
fprintf(fid,"X1_1prime_ L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends MinusiSigma1\n");
case '+igen3sigma_2_SOI'
fprintf(fid,"* +igen3sigma_2_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiGen3Sigma2 L_n1_prime L_n2_prime L_n3_prime R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1 R_n2_primeVarR0  \n");
fprintf(fid,"R2 L_n2 R_n1_prime VarR0  \n");
fprintf(fid,"R3 L_n3 R_n3_prime VarR0  \n");
fprintf(fid,"R4 L_n1 R_n1_prime VarR0_2\n");
fprintf(fid,"R5 L_n2 R_n3_prime VarR0_2\n");
fprintf(fid,"R6 L_n3 R_n2_prime VarR0_2\n");
fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,"R1_ L_n1_prime R_n2  VarR0  \n");
fprintf(fid,"R2_ L_n2_prime R_n1  VarR0  \n");
fprintf(fid,"R3_ L_n3_prime R_n3  VarR0  \n");
fprintf(fid,"R4_ L_n1_prime R_n3  VarR0_2\n");
fprintf(fid,"R5_ L_n2_prime R_n2  VarR0_2\n");
fprintf(fid,"R6_ L_n3_prime R_n1  VarR0_2\n");
fprintf(fid,"X1_1prime_ R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends PlusiGen3Sigma2\n");
case '-igen3sigma_2_SOI'
fprintf(fid,"* -igen3sigma_2_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiGen3Sigma2 L_n1_prime L_n2_prime L_n3_prime R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"R1 L_n1_prime R_n2  VarR0  \n");
fprintf(fid,"R2 L_n2_prime R_n1  VarR0  \n");
fprintf(fid,"R3 L_n3_prime R_n3 VarR0  \n");
fprintf(fid,"R4 L_n1_prime R_n3  VarR0_2\n");
fprintf(fid,"R5 L_n2_prime R_n2 VarR0_2\n");
fprintf(fid,"R6 L_n3_prime R_n1  VarR0_2\n");
fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,"R1_ L_n1 R_n2_primeVarR0  \n");
fprintf(fid,"R2_ L_n2 R_n1_prime VarR0  \n");
fprintf(fid,"R3_ L_n3 R_n3_prime VarR0  \n");
fprintf(fid,"R4_ L_n1 R_n1_prime VarR0_2\n");
fprintf(fid,"R5_ L_n2 R_n3_prime VarR0_2\n");
fprintf(fid,"R6_ L_n3 R_n2_prime VarR0_2\n");
fprintf(fid,"X1_1prime_ L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends MinusiGen3Sigma2\n");
case '+igen3sigma_3_SOI'
fprintf(fid,"* +igen3sigma_3_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt PlusiGen3Sigma3 L_n1_prime L_n2_prime L_n3_prime R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1 R_n1_prime VarC0 \n");
fprintf(fid,"C2 L_n2 R_n2_prime VarC0 \n");
fprintf(fid,"C3 L_n3 R_n3_prime VarC0 \n");
fprintf(fid,"C4 L_n1 R_n2_prime Var2C0\n");
fprintf(fid,"C5 L_n2 R_n3_prime Var2C0\n");
fprintf(fid,"C6 L_n3 R_n1_prime Var2C0\n");
fprintf(fid,"X1_1prime L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,"C1_ L_n1_prime R_n1 VarC0 \n");
fprintf(fid,"C2_ L_n2_prime R_n2 VarC0 \n");
fprintf(fid,"C3_ L_n3_prime R_n3 VarC0 \n");
fprintf(fid,"C4_ L_n1_prime R_n3 Var2C0\n");
fprintf(fid,"C5_ L_n2_prime R_n1 Var2C0\n");
fprintf(fid,"C6_ L_n3_prime R_n2 Var2C0\n");
fprintf(fid,"X1_1prime_ R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends PlusiGen3Sigma3\n");
case '-igen3sigma_3_SOI'
fprintf(fid,"* -igen3sigma_3_SOI \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt MinusiGen3Sigma3 L_n1 L_n2 L_n3 R_n1_prime R_n2_prime R_n3_prime TOGND InitV=0V " + ...
"VarC0=2.7
fprintf(fid,"C1 L_n1_prime R_n1 VarC0 \n");
fprintf(fid,"C2 L_n2_prime R_n2 VarC0 \n");
fprintf(fid,"C3 L_n3_prime R_n3 VarC0 \n");
fprintf(fid,"C4 L_n1_prime R_n3 Var2C0\n");
fprintf(fid,"C5 L_n2_prime R_n1 Var2C0\n");
fprintf(fid,"C6 L_n3_prime R_n2 Var2C0\n");
fprintf(fid,"X1_1prime R_n1_prime R_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime R_n2_prime R_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime R_n3_prime R_n3 TOGND VoltageFollower \n");
fprintf(fid,"C1_ L_n1 R_n1_prime VarC0 \n");
fprintf(fid,"C2_ L_n2 R_n2_prime VarC0 \n");
fprintf(fid,"C3_ L_n3 R_n3_prime VarC0 \n");
fprintf(fid,"C4_ L_n1 R_n2_prime Var2C0\n");
fprintf(fid,"C5_ L_n2 R_n3_prime Var2C0\n");
fprintf(fid,"C6_ L_n3 R_n1_prime Var2C0\n");
fprintf(fid,"X1_1prime_ L_n1_prime L_n1 TOGND VoltageFollower \n");
fprintf(fid,"X2_2prime_ L_n2_prime L_n2 TOGND VoltageFollower \n");
fprintf(fid,"X3_3prime_ L_n3_prime L_n3 TOGND VoltageFollower \n");
fprintf(fid,".ends MinusiGen3Sigma3\n");
end
fprintf(fid,"* ---------------\n");
fprintf(fid,"\n");
end
