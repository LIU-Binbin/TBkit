function H_hr = supercell(H_hr,Ns,filename,Rm,sites,Atom_name,Atom_num,findir)
arguments
    H_hr HR;
    Ns = eye(3);
    filename ='POSCAR_super';
    Rm = H_hr.Rm
    sites = H_hr.sites;
    Atom_name  = H_hr.Atom_name;
    Atom_num = H_hr.Atom_num;
    findir = [0,0,0];
end
format long;
V_Ns=dot(Ns(1,:),cross(Ns(2,:),Ns(3,:)));
 Atom_num_s=round(Atom_num*V_Ns);
[Rm_s, sites_s] = supercell(Ns, Rm, sites, Atom_name, Atom_num, findir, filename);

% ------------------------------------------------------------------------------
H_hr.Atom_num = Atom_num_s;
H_hr.Atom_name = Atom_name;
H_hr.sites=sites_s;
H_hr.Rm=Rm_s;
% H_hr.POSCAR_gen(filename);
end
