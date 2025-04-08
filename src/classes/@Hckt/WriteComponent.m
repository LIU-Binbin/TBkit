function WriteComponent(Component,fid,magnitude,ErrorMode,Error,OPAMP)
arguments
Component char;
fid
magnitude;
ErrorMode = 0;
Error = 0;
OPAMP = '';
end
fprintf(fid,"* ---------------\n");
switch Component
case 'E_A_LC'
fprintf(fid,"* onsite term\n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt E_A TONET TOGND  VarC0=100%s VarL0=1u InitV=0V R_L=1u\n",magnitude);
if ErrorMode
fprintf(fid,"XC0 TONET TOGND Co C_o=GAUSS(VarC0,%7.5f,3) IC=InitV\n",Error);
fprintf(fid,"XL0 TONET TOGND Lo L_o=GAUSS(VarL0,%7.5f,3) R_L=R_L\n",Error);
fprintf(fid,".ends E_A\n");
fprintf(fid,".SubCkt Lo n+ n- L_o=1u R_L=1u\n");
fprintf(fid,"Lp n+ n- L_o R = R_L\n");
fprintf(fid,".ends Lo\n");
else
fprintf(fid,"C0 TONET TOGND VarC0 IC=InitV\n");
fprintf(fid,"L0 TONET TOGND VarL0 R=R_L\n");
fprintf(fid,".ends E_A\n");
end
case 'hopping_C'
fprintf(fid,"* C \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt C n+ n- TOGND  C_hopping=100%s \n",magnitude);
if ErrorMode
fprintf(fid,"Xp n+ n- Co C_o=GAUSS(C_hopping,%7.5f,3)  \n",Error);
fprintf(fid,".ends C\n");
fprintf(fid,".SubCkt Co n+ n- C_o=100p\n");
fprintf(fid,"Cp n+ n- C_o\n");
fprintf(fid,".ends Co\n");
else
fprintf(fid,"Cp n+ n- C_hopping \n");
fprintf(fid,".ends C\n");
end
case 'hopping_minusC_realOPAMP'
fprintf(fid,"* minusC \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt minusC n+ n- TOGND C_hopping = 100%s R_0 = 1 R_L=1u\n",magnitude);
fprintf(fid,"V1 V+ V- 0 DC=3V\n");
fprintf(fid,"X_opamp1 n+ Op1-  V+ V- Op1out %s \n",OPAMP);
fprintf(fid,"X_opamp2 n+ Op2-  V+ V- Op2-  %s \n",OPAMP);
fprintf(fid,"X_opamp3 n- Op3-  V+ V- Op3out %s \n",OPAMP);
fprintf(fid,"X_opamp4 n- Op4-  V+ V- Op4- %s \n",OPAMP);
if ErrorMode
fprintf(fid,"XC_Op1 n+ Op1out Co  C_o=GAUSS(C_hopping,%7.5f,3)\n",Error);
fprintf(fid,"XC_Op3 Op3out n- Co  C_o=GAUSS(C_hopping,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op1 Op1out Op1- Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op3 Op3- Op3out Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op4 Op1- Op4-   Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op2 Op2- Op3-   Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,".ends minusC\n");
fprintf(fid,".SubCkt Ro n+ n- R_o=1000\n");
fprintf(fid,"Rp n+ n- R_o\n");
fprintf(fid,".ends Ro\n");
else
fprintf(fid,"C_Op1 n+ Op1out C_hopping\n");
fprintf(fid,"C_Op3 Op3out n- C_hopping\n");
fprintf(fid,"R_Op1 Op1out Op1- R_0\n");
fprintf(fid,"R_Op3 Op3- Op3out R_0\n");
fprintf(fid,"R_Op4 Op1- Op4- R_0\n");
fprintf(fid,"R_Op2 Op2- Op3- R_0\n");
fprintf(fid,".ends minusC\n");
end
case 'hopping_minusC'
fprintf(fid,"* minusC \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt minusC n+ n- TOGND C_hopping = 100%s R_0 = 1 R_L=1u\n",magnitude);
fprintf(fid,"E_opamp1 Op1out TOGND   n+ Op1-  1E6 max=+100 min=-100\n");
fprintf(fid,"E_opamp2 Op2-   TOGND   n+ Op2-  1E6 max=+100 min=-100\n");
fprintf(fid,"E_opamp3 Op3out TOGND   n- Op3-  1E6 max=+100 min=-100\n");
fprintf(fid,"E_opamp4 Op4-   TOGND   n- Op4-  1E6 max=+100 min=-100\n");
if ErrorMode
fprintf(fid,"XC_Op1 n+ Op1out Co  C_o=GAUSS(C_hopping,%7.5f,3)\n",Error);
fprintf(fid,"XC_Op3 Op3out n- Co  C_o=GAUSS(C_hopping,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op1 Op1out Op1- Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op3 Op3- Op3out Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op4 Op1- Op4-   Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,"XR_Op2 Op2- Op3-   Ro R_o=GAUSS(R_0,%7.5f,3)\n",Error);
fprintf(fid,".ends minusC\n");
fprintf(fid,".SubCkt Ro n+ n- R_o=1000\n");
fprintf(fid,"Rp n+ n- R_o\n");
fprintf(fid,".ends Ro\n");
else
fprintf(fid,"C_Op1 n+ Op1out C_hopping\n");
fprintf(fid,"C_Op3 Op3out n- C_hopping\n");
fprintf(fid,".ends minusC\n");
end
case 'hopping_minusC_DM'
fprintf(fid,"* minusC 2 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt minusC_DM n+ n- TOGND C_hopping = 100%s C_h = 'SQRT(C_hopping)' L_DM = 1u C_DM = 1\n",magnitude);
fprintf(fid,"Xonsite n1 TOGND E_A VarC0 =  C_DM  VarL0 = L_DM\n");
fprintf(fid,"C1 n+ n1 C_h\n");
fprintf(fid,"C2 n1 n- C_h\n");
fprintf(fid,".ends minusC_DM\n");
case 'E_A_RC'
fprintf(fid,"* onsite term\n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt E_A TONET TOGND  VarR0=1 VarC0=100%s InitV=0V \n",magnitude);
fprintf(fid,"R0 TONET TOGND VarR0\n");
fprintf(fid,"C0 TONET TOGND VarC0 IC=InitV\n");
fprintf(fid,".ends E_A\n");
case 'hopping_L'
fprintf(fid,"* L \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt L n+ n- TOGND  L_hopping=1%s R_L=1u\n",magnitude);
fprintf(fid,"Lp n+ n- L_hopping R=R_L\n");
fprintf(fid,".ends L\n");
case 'hopping_AL'
fprintf(fid,"* hopping_AL \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt AL n+ n- TOGND L_hopping2 = 2%s R_L=1u\n",magnitude);
fprintf(fid,"E_opamp1 Op1-   TOGND   n+ Op1-  level=1\n");
fprintf(fid,"Lp2 Op1- n- L_hopping2 R=R_L\n");
fprintf(fid,".ends AL\n");
case 'hopping_LA'
fprintf(fid,"* hopping_LA \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt LA n+ n- TOGND L_hopping2 = 2%s R_L=1u\n",magnitude);
fprintf(fid,"E_opamp1 Op1-   TOGND   n- Op1-  level=1 \n");
fprintf(fid,"Lp2 Op1- n+ L_hopping2 R=R_L\n");
fprintf(fid,".ends LA\n");
case 'Port_L3'
fprintf(fid,"* Lport3 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt Lport3 TONET TOGND  L_hopping=1%s L_hopping2=2%s L_hopping3=3%s R_L=1u \n",magnitude,magnitude,magnitude);
fprintf(fid,"Lp1 TONET TOGND L_hopping R=R_L\n");
fprintf(fid,"Lp2 TONET TOGND L_hopping2 R=R_L\n");
fprintf(fid,"Lp3 TONET TOGND L_hopping3 R=R_L\n");
fprintf(fid,".ends Lport3\n");
case 'Port_L1'
fprintf(fid,"* Lport1 \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt Lport1 TONET TOGND  L_hopping=1%s R_L=1u\n",magnitude);
fprintf(fid,"Lp TONET TOGND L_hopping R=R_L \n");
fprintf(fid,".ends Lport1\n");
end
fprintf(fid,"* ---------------\n");
fprintf(fid,"\n");
end
