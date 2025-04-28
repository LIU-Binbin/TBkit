function [orbital_out,Rm_s_fin] = AddVacuumLayer(orbital_init,POSCAR_file,fin_dir_list,options)
            arguments
                orbital_init
                POSCAR_file = 'POSCAR';
                fin_dir_list = [1 1 0];
                options.fast = true;
                options.vacuum_length = 10;
            end
            orbital_out = orbital_init;
            fin_orb = orbital_out;
            Ns = [1 0 0;0 1 0;0 0 1];
            if ~options.fast
                %disp(fin_dir_list);
                [Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=vasplib.POSCAR_read(POSCAR_file);
                % gen POSCAR
                H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
            else
                switch class(POSCAR_file)
                    case {'string','char'}
                        [Rm_tmp,~,~,~]=vasplib.POSCAR_read(POSCAR_file);
                    case 'double'
                        Rm_tmp = POSCAR_file;
                end
            end
            % rebuild fin_orb
            Rm_tmp = Ns*Rm_tmp;
            Rmlength1 = norm (Rm_tmp(1,:));
            Rmlength2 = norm (Rm_tmp(2,:));
            Rmlength3 = norm (Rm_tmp(3,:));
            Rm_s_fin_add = [options.vacuum_length*Rm_tmp(1,:)*fin_dir_list(1)/Rmlength1;...
                options.vacuum_length*Rm_tmp(2,:)*fin_dir_list(2)/Rmlength2;...
                options.vacuum_length*Rm_tmp(3,:)*fin_dir_list(3)/Rmlength3];
            Rm_s_fin = Rm_tmp + Rm_s_fin_add ;
            Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
            Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
            [nfinorb,~ ]= size(fin_orb);
            for  i = 1:nfinorb
                Rr_orb = fin_orb(i,:)*Rm_tmp;
                Rr_s_fin = Rr_orb + Rr_s_fin_add;
                Rc_s_fin = Rr_s_fin / Rm_s_fin;
                fin_orb(i,:) = Rc_s_fin ;
            end
            orbital_out = fin_orb;% change coordinate along finite direction ; fractional
        end
