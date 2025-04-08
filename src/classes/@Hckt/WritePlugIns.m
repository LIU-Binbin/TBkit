function WritePlugIns(PlugIns,fid,magnitude)
arguments
PlugIns char;
fid;
magnitude = 'p';
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
switch PlugIns
case 'Na-(Mb+Lc)'
fprintf(fid,"* Na-(Mb+Lc) \n");
fprintf(fid,"*\n");
fprintf(fid,"*Rb = Rf/M; Rc = Rf/L Ra = (2+M+L-N)Rf RNf = N Rf \n");
fprintf(fid,".SubCkt AdderSubtractor va vb vc vo TOGND " + ...
"VarRf=806
Rmagnitude,Rmagnitude,Rmagnitude,Rmagnitude,Rmagnitude);
fprintf(fid,"Ra va v_plus VarRa \n");
fprintf(fid,"Rb vb v_minus VarRb \n");
fprintf(fid,"Rc vc v_minus VarRc \n");
fprintf(fid,"Rf1 TOGND v_minus VarRf \n");
fprintf(fid,"Rf2 vo v_minus VarRf \n");
fprintf(fid,"RNf TOGND v_plus VarRNf \n");
fprintf(fid,"E_opamp vo TOGND v_plus v_minus  level=1\n");
fprintf(fid,".ends AdderSubtractor\n");
case 'VoltageFollower'
fprintf(fid,"* VoltageFollower \n");
fprintf(fid,"*\n");
fprintf(fid,".SubCkt VoltageFollower ui uo TOGND\n");
fprintf(fid,"E_opamp uo TOGND ui uo  level=1\n");
fprintf(fid,".ends VoltageFollower\n");
case 'InvertingOpAmp'
fprintf(fid,"*InvertingOpAmp (M1) \n");
fprintf(fid,"*\n");
fprintf(fid,"* va vb TOGND \n");
fprintf(fid,".SubCkt InvertingOpAmp va vb TOGND " + ...
"VarRf=7.87k
fprintf(fid,"R1 vb v_minus Var2Rf \n");
fprintf(fid,"R2 va v_minus VarRf\n");
fprintf(fid,"E_opamp vb TOGND TOGND v_minus level=1\n");
fprintf(fid,".ends InvertingOpAmp\n");
case 'IntegratorOpAmp'
fprintf(fid,"*IntegratorOpAmp (M2) \n");
fprintf(fid,"*\n");
fprintf(fid,"* va vo TOGND \n");
fprintf(fid,".SubCkt IntegratorOpAmp va vo TOGND " + ...
" VarCm=47
%fprintf(fid,"C1 vo vb Var3Cm \n");
fprintf(fid,"C2 vo v_minus VarCm\n");
fprintf(fid,"R1 vo v_minus VarR \n");
fprintf(fid,"R2 va v_minus VarRm\n");
fprintf(fid,"R3 TOGND v_plus VarRm\n");
fprintf(fid,"E_opamp vo TOGND v_plus v_minus level=1\n");
fprintf(fid,".ends IntegratorOpAmp\n");
end
end
