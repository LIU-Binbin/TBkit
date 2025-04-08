function [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_gen_3D(H_hr,kmesh,k3d,options)
arguments
H_hr HR
kmesh = [100,100];
k3d = [-0.5, -0.5, 0;...
1.0 , 0,   0;...
0   , 1.0,   0;...
];
options.output ='raw_data';
options.LWAVE = false;
options.cartisian = true;
options.fin_dir = 3;
options.ProjectionMethod = 'slab';
options.ProjectionStruct = struct('discrimination',0.1,'center',[0.5,0.5,0.5],'orientation',2,'sign',false);
options.WEIGHTCAR = false;
options.norb = -1;
options.fermi = 0;
end
options2 = rmfield(options,{'output','cartisian','fin_dir'});
optionsCell = namedargs2cell(options2);
[H_hr.klist_frac,klist1,klist2] = H_hr.kmesh3D(kmesh,k3d,'fermiarc');
if options.WEIGHTCAR
[EIGENCAR_3D,WAVECAR_3D,WEIGHTCAR_3D] = H_hr.EIGENCAR_gen(...
optionsCell{:});
else
[EIGENCAR_3D,WAVECAR_3D] = H_hr.EIGENCAR_gen('LWAVE',options.LWAVE);
WEIGHTCAR_3D = [];
end
if options.LWAVE
else
WAVECAR_3D = [];
end
if options.cartisian
options.output = 'refined';
klist = H_hr.klist_frac*H_hr.Gk;
end
if options.LWAVE
WEIGHTCAR_3D = WAVECAR_3D;
end
if strcmp(options.output ,'refined')
switch options.fin_dir
case 1
klist1 =klist(:,2);
klist2 =klist(:,3);
case 2
klist1 =klist(:,1);
klist2 =klist(:,3);
case 3
klist1 =klist(:,1);
klist2 =klist(:,2);
end
klist1 = reshape(klist1,kmesh);
klist2 = reshape(klist2,kmesh);
EIGENCAR_3D = reshape(EIGENCAR_3D.',kmesh(1),kmesh(2),[]);
if options.WEIGHTCAR
WEIGHTCAR_3D = reshape(WEIGHTCAR_3D.',kmesh(1),kmesh(2),[]);
end
end
end
