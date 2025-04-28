        function [Rm,sites,Atom_name,Atom_num] = POSCAR_gen(TBkitobj,filename,Rm,sites,Atom_name,Atom_num)
            %% note: POSCAR is Direct mode
            % input : Stucture information
            %                             a_crystal_constance
            %                             a1,a2,a3,Rm=[a1;a2;a3]
            %                             atom_name;atom_num
            %                             sites
            %          element_information:
            %                             sites
            % Output : POSCAR:
            %                             a_crystal_constance
            %                             a1,a2,a3,Rm=[a1;a2;a3]
            %                             atom_name;atom_num
            %                             coordinates pattern

            % note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
            % Rm
            %a1 a1x a1y a1z
            %a2 a2x a2y a2z
            %a3 a3x a3y a3z
            % unit A?e-6
            % usage : POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);
            arguments
                TBkitobj
                filename = 'POSCAR';
                Rm = TBkitobj.Rm;
                sites = TBkitobj.sites;
                Atom_name = TBkitobj.Atom_name;
                Atom_num = TBkitobj.Atom_num;
            end
            title = "POSCAR Gen by vasplib";
            %% write POSCAR
            % Initialize variables
            %filename = "POSCAR_"+num2str(term2)+".vasp";
            fileID = fopen(filename,'w');
            %% Rm
            a_crystal_constance=1;
            %%
            fprintf(fileID,"%s\n",title);
            fprintf(fileID,"%d\n",a_crystal_constance);
            %fprintf(fileID,"  ",Rm(i,j));
            for i=1:3
                for j=1:3
                    fprintf(fileID,"  %f",Rm(i,j));
                end
                fprintf(fileID,"\n");
            end
            for i=1:length(Atom_name)
                fprintf(fileID,"%s ",Atom_name(i));
            end
            fprintf(fileID,"\n  ");
            for i=1:length(Atom_num)
                fprintf(fileID,"%d ",Atom_num(i));
            end
            fprintf(fileID,"\n");
            fprintf(fileID,"Direct\n  ");
            % sites
            [~,sites_num]=size(sites);
            for i=1:sites_num
                fprintf(fileID,"%f  ",mod(sites(i).rc1,1));
                fprintf(fileID,"%f  ",mod(sites(i).rc2,1));
                fprintf(fileID,"%f  ",mod(sites(i).rc3,1));
                fprintf(fileID,"\n  ");
            end
            fclose(fileID);
        end

