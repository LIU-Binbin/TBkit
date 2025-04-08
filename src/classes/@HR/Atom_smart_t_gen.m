function Atom_smart_t = Atom_smart_t_gen(site1,site2)
Rc1 = [site1.rc1,site1.rc2,site1.rc3];
Rc2 = [site2.rc1,site2.rc2,site2.rc3];
Atom_smart_t.R_fractional_from = Rc1 ;
Atom_smart_t.R_fractional_to   = Rc2;
Atom_smart_t.R_fractional_diff = -(Rc1 - Rc2);
Atom_smart_t.seq_from = site1.seq;
Atom_smart_t.seq_to = site2.seq;
Atom_smart_t.l_name_from = site1.orb  ;
Atom_smart_t.l_name_to   = site2.orb  ;
Atom_smart_t.orb_sym_from = site1.orb_sym;
Atom_smart_t.orb_sym_to = site2.orb_sym;
Atom_smart_t.handyname = strcat(site1.name,' -> ',site2.name);
end
