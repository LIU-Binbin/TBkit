function EIGENCAR = ProjectDOSCAR(DOSCAR,options)
arguments
DOSCAR;
options.POSCAR = 'POSCAR';
options.KPOINTS = 'KPOINTS';
options.dir_seq = [1,2,3];
options.klist_frac = [];
options.method{mustBeMember(options.method,{'linear','nearest','pchip','cubic','makima','spline'})} = 'nearest';
end
if exist(options.KPOINTS,'file') && exist(options.POSCAR,'file')
Rm=POSCAR_read(options.POSCAR);
[~,~,klist_frac,~,~]=kpathgen3D(Rm,options.KPOINTS);
else
if isempty(options.klist_frac)
if length(size(DOSCAR)) >2
error('No KPOINTS POSCAR');
else
klist_frac = zeros(size(DOSCAR,2),1);
end
else
klist_frac = options.klist_frac;
end
end
if isa(DOSCAR,'cell')
NBANDS = length(DOSCAR);
sizemesh = size(DOSCAR{1});
Dim = length(sizemesh);
cellmode = true;
else
sizemesh = size(DOSCAR);
Dim = length(sizemesh)-1;
if Dim == 1
NBANDS = sizemesh(1);
else
NBANDS = sizemesh(Dim+1);
end
cellmode = false;
end
if length(options.dir_seq) ~= Dim
options.dir_seq = 1:Dim;
end
klist_frac = mod(klist_frac,1);
nklist = size(klist_frac,1);
if Dim > 3 || cellmode
EIGENCAR = zeros(NBANDS,nklist);
klistd_frac{Dim} = linspace(0,1,sizemesh(Dim));
KNgrid{Dim} = [];
klist_fracCell{Dim} = [];
for d = 1:Dim
klistd_frac{d} = linspace(0,1,sizemesh(d));
klist_fracCell{d} = klist_frac(:,options.dir_seq(d));
end
[KNgrid{:}] = ndgrid(klistd_frac{:});
for i = 1:NBANDS
V = DOSCAR{i};
EIGENCAR(i,:) = interpn(KNgrid{:},V,klist_fracCell{:},options.method).';
end
elseif Dim == 1
EIGENCAR = zeros(NBANDS,nklist);
klist1_s = linspace(0,1,sizemesh(2));
TargetMat = normalize(DOSCAR,'scale');
TargetMatMax =  max(DOSCAR,[],'all');
TargetMatMin =  min(DOSCAR,[],'all');
TargetMatPlus = DOSCAR > 0;
TargetMatMinus = DOSCAR < 0;
TargetMatABSMax = max(abs(TargetMatMax),abs(TargetMatMin));
TargetMatPlusFactor = TargetMatMax/TargetMatABSMax;
TargetMatMinusFactor = TargetMatMin/TargetMatABSMax;
if ~sum(TargetMatPlus,'all')
TargetMat(TargetMatPlus) = normalize(TargetMat(TargetMatPlus),'range',[0,1]);
end
if ~sum(TargetMatMinus,'all')
TargetMat(TargetMatMinus) = normalize(TargetMat(TargetMatMinus),'range',[-1,0]);
end
for i = 1:sizemesh(1)
EIGENCAR(i,:) = interp1(klist1_s,TargetMat(i,:),klist_frac(:,options.dir_seq(1)),options.method);
end
elseif Dim == 2
EIGENCAR = zeros(NBANDS,nklist);
klist1_s = linspace(0,1,sizemesh(1));
klist2_s = linspace(0,1,sizemesh(2));
[K1grid,K2grid] = meshgrid(klist1_s,klist2_s);
for i = 1:NBANDS
EIGENCAR(i,:) = interp2(K1grid,K2grid,DOSCAR(:,:,i),klist_frac(:,options.dir_seq(1)),klist_frac(:,options.dir_seq(2)),options.method);
end
elseif Dim == 3
EIGENCAR = zeros(NBANDS,nklist);
klistd_frac{Dim} = linspace(0,1,sizemesh(Dim));
KNgrid{Dim} = [];
klist_fracCell{Dim} = [];
for d = 1:Dim
klistd_frac{d} = linspace(0,1,sizemesh(d));
klist_fracCell{d} = klist_frac(:,options.dir_seq(d));
end
[KNgrid{:}] = ndgrid(klistd_frac{:});
for i = 1:NBANDS
V = DOSCAR(:,:,:,i);
EIGENCAR(i,:) = interpn(KNgrid{:},V,klist_fracCell{:},options.method).';
end
else
end
end
