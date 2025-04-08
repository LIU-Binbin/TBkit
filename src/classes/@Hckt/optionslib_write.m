function optionslib_write(HcktObj,fid,options3,options5)
arguments
HcktObj
fid
options3.libmode = false;
options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
end
fprintf(fid,"* -------- Initial --------\n");
for i = 1:size(HcktObj.IC,1)
fprintf(fid,'* ');
fprintf(fid,(HcktObj.IC(i,:)));fprintf(fid,"\n");
end
fprintf(fid,"* -------- OPTIONS --------\n");
switch  HcktObj.magnitude
case 'p'
if strcmp(HcktObj.TRAN,"")
HcktObj.TRAN =".tran 1ns 10us";
end
if strcmp(HcktObj.AC,"")
HcktObj.AC = ".AC LIN 1000 0.5e07 2e07";
end
DEFAULT_PARM = "
case 'n'
if strcmp(HcktObj.TRAN,"")
HcktObj.TRAN =".tran 50ns 500us";
end
if strcmp(HcktObj.AC,"")
HcktObj.AC =".AC LIN 1000 0.1e06 5e06";
end
DEFAULT_PARM = "
case 'u'
if strcmp(HcktObj.TRAN,"")
HcktObj.TRAN =".tran 1us 10ms";
end
if strcmp(HcktObj.AC,"")
HcktObj.AC = ".AC LIN 1000 0.5e04 2e04";
end
DEFAULT_PARM = "
end
switch options5.analysis
case 'ac'
OPTcase = [HcktObj.AC;HcktObj.Options;DEFAULT_PARM];
case 'tran'
OPTcase = [HcktObj.TRAN;HcktObj.Options;DEFAULT_PARM];
otherwise
OPTcase = [HcktObj.TRAN;HcktObj.Options;DEFAULT_PARM];
end
for i = 1:size(OPTcase,1)
fprintf(fid,(OPTcase(i,:)));fprintf(fid,"\n");
end
fprintf(fid,"* -------- -------- --------\n");
fprintf(fid,"* -------- lib --------\n");
for i = 1:size(HcktObj.Lib,1)
fprintf(fid,(HcktObj.Lib(i,:)));fprintf(fid,"\n");
end
fprintf(fid,"\n");
fprintf(fid,"* -------- -------- --------\n");
fprintf(fid,"* -------- Subckt for Homecell --------\n");
fprintf(fid,"*\n");
fprintf(fid,HcktObj.HomeCell.OutputStr);
fprintf(fid,"\n");
fprintf(fid,"*\n");
fprintf(fid,"* -------- Subckt for Hopping --------\n");
fprintf(fid,"*\n");
if options3.libmode
for i = 2: HcktObj.NRPTS
if HcktObj.ScktL(i).type == 'X'
fprintf(fid,"* %s\n",mat2str(HcktObj.SktL(i).name));
fprintf(fid,HcktObj.ScktL(i).OutputStr);
fprintf(fid,"\n");
fprintf(fid,"* --------");
end
end
end
fprintf(fid,"*\n");
switch options5.analysis
case 'tran'
fprintf(fid,"* -------- Ipulse --------\n");
fprintf(fid,"%s\n",HcktObj.Ipulse);
case 'ac'
fprintf(fid,"* -------- V ac --------\n");
fprintf(fid,"%s\n",HcktObj.Vac);
otherwise
end
end
