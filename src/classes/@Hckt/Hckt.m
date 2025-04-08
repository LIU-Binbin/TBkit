classdef Hckt < matlab.mixin.CustomDisplay
properties
title = "Created by Hckt "+date;
Nports ;
MeshNnodes ;
vectorAll = [0,0];
end
properties
HnodeL  logical= [];
ScktL Subckt = Subckt([]);
vectorL = [1];
PortInCell = {[]};
PortOutCell = {[]};
DescriptionL = "";
end
properties
NetlistMesh;
fin_dir = [];
mesh = [];
Port2VectorDist ;
Port2PortDist ;
Port2PortDistForVectorise ;
Port2PortDist_BinBin ;
AddtionInformation ='';
end
properties
magnitude = 'p';
Lib = [".inc 'Hckt.lib' ";];
IC = ["";];
Ipulse = ["";];
Vac =["";];
TRAN = "";
AC = "";
Options = [...
".option post=2 probe";...
"*.option parhier=global";...
".option accurate=1";...
".option MTTHRESH=2";...
".option FAST";...
"*.option INTERP";]
end
properties
VectorizeScktL Subckt = Subckt([]) ;
end
properties
pinned = false;
end
properties(Dependent)
NRPTS;
HomeCell Subckt  ;
dim;
Components;
innerNetlist;
NetlistStr;
Netlist;
pinnedNetlist;
OPT;
end
properties
pinnedNetlistStr;
end
properties(Transient =true)
end
methods (Access = protected)
 propgrp = getPropertyGroups(~)
end
methods
 Hcktobj = Hckt(ScktL,vectorAll,propArgs,options)
end
methods
 Hcktobj = set_home(Hcktobj,Subcktobj,PortInL,PortOutL,DescriptionL)
 Hcktobj = set_hop(Hcktobj,vector,Subcktobj,PortInL,PortOutL,DescriptionL,options)
 Hcktobj = add_empty_one(Hcktobj,vector)
end
methods
 OPT = get.OPT(HcktObj)
 HomeCell = get.HomeCell(HcktObj)
 Components = get.Components(HcktObj)
 NRPTS = get.NRPTS(HcktObj)
 dim = get.dim(HcktObj)
 innerNetlist = get.innerNetlist(HcktObj)
 NetlistStr = get.NetlistStr(HcktObj)
 pinnedNetlist = get.pinnedNetlist(HcktObj)
 Netlist = get.Netlist(HcktObj)
end
methods(Static)
 [simulation_result]=read_hspice_ac(filename,options)
 [simulation_result]=read_hspice_tr(filename,options)
 [simulation_result]=read_hspice_tr_sw_ac(filename)
end
methods(Static)
 simulation_result = save_signal_names(num_var,data_start_ind,content_str,file_extension)
end
methods(Static)
 [VectorList,ObservationsMat,Vstruct,SpectrumX,TimeL] = extractObservations(simulation_result,Observations,options)
 EIGENCAR = ProjectDOSCAR(DOSCAR,options)
 [DOSCAR_2D,klist1,OmgL] = CollectVstruct1D(ObservationsMat,VectorList,OmegaCut,SpectrumX)
 [DOSCAR_3D,klist1,klist2,OmgL] = CollectVstruct2D(VectorList,ObservationsMat,OmegaCut,SpectrumX)
 [DOSCAR_4D,klist1,klist2,klist3,OmgL] = CollectVstruct3D(VectorList,ObservationsMat,OmegaCut,SpectrumX,options)
 [DOSCAR_ND,klistCell,OmgL] = CollectVstruct(VectorList,ObservationsMat,OmegaCut,SpectrumX,options)
end
methods(Static)
 [fig,ax] = Frequencyplot(ObservationsMat,SpectrumL,OmegaCut,options)
 [ax,WaveFunc] = waveplot(ObservationsMat,VectorList,SpectrumL,orbL,options_select,options)
 [ax] = bandplot(EIGNECAR,OmegaCut,SpectrumL,options,optionskpath)
end
methods
 HcktObj = vectorize(HcktObj,fin_dir)
 [HcktObj,Basename] = hspice_gen(HcktObj,filename,ops,options,options2,options3,options4,options5)
 optionslib_write(HcktObj,fid,options3,options5)
 HcktObj = hspice_gen_vectorized(HcktObj,filename,options,options2,options3)
 HcktObj = hspice_gen_general(HcktObj,filename,options,options2,options3,options5)
 [ScktnameL] = hoppinglist_gen(HcktObj,fid,seq,ijkL,ScktnameL,options)
 [HcktObj,ScktnameL,AllPortStrL,ijkL] = netlist_gen(HcktObj,fid,options)
end
methods(Static)
 WriteVendorComponentLibrary
 WritePlugIns(PlugIns,fid,magnitude)
 WriteModules(Modules,fid,magnitude,Lprefix)
 WriteComponent(Component,fid,magnitude,ErrorMode,Error,OPAMP)
 Genlib(filename,options)
 Portlist_gen(fid,AllPortStrL,options,options2)
end
methods
 HcktObj = autohermi(HcktObj)
 HcktObj = half(HcktObj,checkcommute)
 HcktObj = reseq(HcktObj,seqL)
 HcktObj = clean(HcktObj)
end
methods(Static)
 [VectorList,ObservationsMat,SpectrumX] = FakeSimulation(H_hr,meshs,OmegaCut,options)
end
methods(Static)
 vectorL = dim2vectorL(dim,maxR)
end
end
