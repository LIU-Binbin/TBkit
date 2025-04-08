function HcktObj = hspice_gen_general(HcktObj,filename,options,options2,options3,options5)
arguments
HcktObj Hckt;
filename = strcat("Hckt_test",".sp");
options.mesh = repmat(10,[HcktObj.dim,1]);
options.fin_dir = zeros(HcktObj.dim,1);
options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
options2.prefix = "+ ";
options2.father = "";
options2.nodelist = [];
options3.fft = false;
options3.libmode = false;
options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
end
optionscell = namedargs2cell(options);
fid = fopen(filename,'w');
fprintf(fid,HcktObj.title);
fprintf(fid,"\n");
HcktObj.optionslib_write(fid,"analysis",options5.analysis );
fprintf(fid,"* -------- Netlist Homecell --------\n");
fprintf(fid,"*\n");
[HcktObj,ScktnameL,AllPortStrL,ijkL] = netlist_gen(HcktObj,fid,optionscell{:});
fprintf(fid,"*\n");
fprintf(fid,"* -------- Netlist Homecell end --------\n");
fprintf(fid,"*\n");
fprintf(fid,"* -------- Netlist Hopping --------\n");
fprintf(fid,"*\n");
for i = 2:HcktObj.Components
fprintf(fid,"* \\ %s|%s <-- %s \\ %s\n",...
mat2str(HcktObj.vectorAll(HcktObj.vectorL(i),:)),...
mat2str(HcktObj.PortOutCell{i,:}),...
mat2str(HcktObj.PortInCell{i,:}),...
HcktObj.ScktL(i).name(1) ...
);
hoppinglist_gen(HcktObj,fid,i,ijkL,ScktnameL,optionscell{:});
end
fprintf(fid,"*\n");
fprintf(fid,"* -------- Netlist Hopping end --------\n");
fprintf(fid,"*\n");
switch options5.analysis
case 'tran'
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
case 'ac'
fprintf(fid,"* -------- Probe List --------\n");
fprintf(fid,"*\n");
fprintf(fid,".PROBE AC\n");
optionsProbe = options2;
optionsProbe.comment = false;
optionsProbe.type = 'VR';
optionsProbeCell = namedargs2cell(optionsProbe);
switch options2.probenode
case 'Allnode'
Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
case 'Selectnode'
Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
end
optionsProbe.type = 'VI';
optionsProbeCell = namedargs2cell(optionsProbe);
switch options2.probenode
case 'Allnode'
Hckt.Portlist_gen(fid,AllPortStrL,optionsProbeCell{:});
case 'Selectnode'
Hckt.Portlist_gen(fid,ScktnameL,optionsProbeCell{:});
end
fprintf(fid,"* -------- Probe List end --------\n");
fprintf(fid,"*\n");
end
fprintf(fid,".END\n");
fclose(fid);
end
