function H_hr = supercell(H_hr,Ns,filename,Rm,sites,Atom_name,Atom_num,findir)
Accuracy=8;
if nargin < 4
Rm = H_hr.Rm;
end
if nargin < 5
sites = H_hr.sites;
end
if nargin < 6
Atom_name = H_hr.Atom_name;
end
if nargin < 7
Atom_num = H_hr.Atom_num;
end
if nargin < 8
findir = [0,0,0];
end
if nargin == 4
Origin_POSCAR = Rm;
[Rm,sites,Atom_name,Atom_num,~]=HR.POSCAR_read(Origin_POSCAR,'vasp');
findir = [0,0,0];
end
if nargin < 3
filename  = 'POSCAR_super';
end
if nargin < 2
Ns = eye(3);
end
format long;
V_Ns=dot(Ns(1,:),cross(Ns(2,:),Ns(3,:)));
if V_Ns < 0
Ns(3,:) = -Ns(3,:);
end
V_Ns = abs(V_Ns);
if V_Ns < 1
error('check Ns');
end
[~,nsites] = size(sites);
for i =1:nsites
sites(i).rc1 = HR.plusrc(sites(i).rc1);
if sites(i).rc1  < 0
disp('bugtest')
end
sites(i).rc2 = HR.plusrc(sites(i).rc2);
if sites(i).rc2  < 0
disp('bugtest')
end
sites(i).rc3 = HR.plusrc(sites(i).rc3);
if sites(i).rc3  < 0
disp('bugtest')
end
end
if V_Ns >= 1
max_R=max(abs(Ns))*3;
sc_cands = [];
for i = -max_R(1):max_R(1)
for j = -max_R(2):max_R(2)
for k = -max_R(3):max_R(3)
sc_cands=[sc_cands;[i,j,k]];
end
end
end
sc_vec=[];
eps_shift=sqrt(2.0)*1.0E-8;
for ivec = 1:length(sc_cands)
tmp_red=HR.to_red_sc(sc_cands(ivec,:),Ns);
inside = 1;
for it = 1:length(tmp_red)
t = tmp_red(it);
if t <= -1.0*eps_shift || t>1.0-eps_shift
inside=0;
end
end
if inside == 1
sc_vec=[sc_vec;sc_cands(ivec,:)];
end
end
[num_sc,~] = size(sc_vec);
if round(round(abs(det(Ns)))) ~= num_sc
error("\n\nSuper-cell generation failed! Wrong number of super-cell vectors found.");
end
countnum=0;
countnum2=0;
for elementseq=1:length(Atom_num)
for n=1:Atom_num(elementseq)
countnum2=countnum2+1;
Rpc1=sites(countnum2).rc1;
Rpc2=sites(countnum2).rc2;
Rpc3=sites(countnum2).rc3;
Rc=[Rpc1 Rpc2 Rpc3];
for icur_sc_vec = 1:num_sc
cur_sc_vec = sc_vec(icur_sc_vec,:);
Rsc = HR.to_red_sc(Rc+cur_sc_vec,Ns);
Rsc = round(Rsc.*10^Accuracy)./10^Accuracy;
rsc1=Rsc(1);
rsc2=Rsc(2);
rsc3=Rsc(3);
if rsc1 ==1
rsc1 = 0;
end
if rsc2 ==1
rsc2 = 0;
end
if rsc3 ==1
rsc3 = 0;
end
countnum=countnum+1;
sites_s(countnum).rc1=rsc1;
sites_s(countnum).rc2=rsc2;
sites_s(countnum).rc3=rsc3;
sites_s(countnum).name=Atom_name(elementseq);
end
end
end
sites=sites_s;
end
Atom_num_s=round(Atom_num*V_Ns);
Rm_s=Ns*Rm;
[~,nsites] = size(sites);
if findir(1) ~=0 || findir(2) ~=0 || findir(3) ~=0
Rm = Ns*Rm;
disp('findir_mode');
Rmlength1 = abs(norm (Rm(1,:)));
Rmlength2 = abs(norm (Rm(2,:)));
Rmlength3 = abs(norm (Rm(3,:)));
Rm_s_fin_add = [10*Rm(1,:)*findir(1)/Rmlength1;...
10*Rm(2,:)*findir(2)/Rmlength2;...
10*Rm(3,:)*findir(3)/Rmlength3];
Rm_s_fin = Rm + Rm_s_fin_add ;
Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
for  i = 1:nsites
Rr = [sites(i).rc1,sites(i).rc2,sites(i).rc3]*Rm;
Rr_s_fin = Rr + Rr_s_fin_add;
Rc_s_fin = Rr_s_fin /Rm_s_fin;
sites_s(i).rc1 = Rc_s_fin(1);
sites_s(i).rc2 = Rc_s_fin(2);
sites_s(i).rc3 = Rc_s_fin(3);
end
Rm_s=Rm_s_fin;
end
H_hr.Atom_num = Atom_num_s;
H_hr.Atom_name = Atom_name;
H_hr.sites=sites_s;
H_hr.Rm=Rm_s;
H_hr.POSCAR_gen(filename);
end
