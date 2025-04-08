function HcktObj = hspice_gen_vectorized(HcktObj,filename,options,options2,options3)
arguments
HcktObj Hckt;
filename = strcat("Hckt_test",".sp");
options.mesh = repmat(10,[HcktObj.dim,1]);
options.fin_dir = zeros(HcktObj.dim,1);
options.BinBin = false;
options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
options2.prefix = "+ ";
options2.father = "";
options2.nodelist = [];
options3.fft = false;
options3.libmode = false;
end
optionscell = namedargs2cell(options);
HcktObj = HcktObj.vectorize(options.fin_dir);
fid = fopen(filename,'w');
fprintf(fid,HcktObj.title);
fprintf(fid,"\n");
HcktObj.optionslib_write(fid,'libmode',options3.libmode);
fprintf(fid,"* -------- Subckt for Vectorize --------\n");
fprintf(fid,"*\n");
fprintf(fid,"*%s\n",HcktObj.VectorizeScktL(1).name);
fprintf(fid,HcktObj.VectorizeScktL(1).OutputStr);
fprintf(fid,"\n");
fprintf(fid,"*\n");
fprintf(fid,"* -------- Subckt for Vectorize (OpenBoundary) --------\n");
fprintf(fid,"*\n");
for i = 2: length(HcktObj.VectorizeScktL)
fprintf(fid,"* %s\n",HcktObj.VectorizeScktL(i).name);
fprintf(fid,HcktObj.VectorizeScktL(i).OutputStr);
fprintf(fid,"\n");
fprintf(fid,"* --------");
end
fprintf(fid,"* -------- -------- --------\n");
fprintf(fid,"* -------- Netlist --------\n");
fprintf(fid,"*\n");
[HcktObj,ScktnameL,AllPortStrL] = HcktObj.netlist_gen(fid,optionscell{:},'convention','II');
fprintf(fid,"*\n");
fprintf(fid,"* -------- Netlist end --------\n");
fprintf(fid,"*\n");
fprintf(fid,"* -------- Probe List --------\n");
fprintf(fid,"*\n");
fprintf(fid,".PROBE TRAN\n");
optionsProbe = options2;
optionsProbe.comment = false;
optionsProbeCell = namedargs2cell(optionsProbe);
switch options2.probenode
case 'Allnode'
Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
case 'Selectnode'
Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
end
fprintf(fid,"* -------- Probe List end --------\n");
fprintf(fid,"*\n");
fprintf(fid,"* -------- FFT List --------\n");
fprintf(fid,"*\n");
optionsFFT = options2;
optionsFFT.comment = ~options3.fft;
optionsFFT.prefix = ".FFT ";
optionsFFTcell = namedargs2cell(optionsFFT);
switch options2.probenode
case 'Allnode'
Hckt.Portlist_gen(fid,AllPortStrL,optionsFFTcell{:});
case 'Selectnode'
Hckt.Portlist_gen(fid,ScktnameL,optionsFFTcell{:});
end
fprintf(fid,"* -------- FFT List end --------\n");
fprintf(fid,"*\n");
fprintf(fid,".END\n");
fclose(fid);
end
