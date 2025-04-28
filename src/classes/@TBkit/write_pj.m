function fid = write_pj(TBkitobj,mode,filename)
            arguments
                TBkitobj ;
                mode = 'wt';
                filename = 'wt.in';
            end
            % function string
            if strcmp(mode,'wt')
                fid = fopen(filename,'a');
                %% wt.in projectors card gen
                maprule3= containers.Map({-1,0,1,2,4,5,6,8},{"","s", "pz px py","dz2 dxz dyz dx2-y2 dxy","s pz px py","s dz2 dxz dyz dx2-y2 dxy","pz px py dz2 dxz dyz dx2-y2 dxy","s pz px py dz2 dxz dyz dx2-y2 dxy"});
                maprule4= containers.Map({-1,0,1,2,3,4,5,6,7,8,9},{0,1,3,5,7,4,6,8,8,9,16});
                % function string
                head_string="PROJECTORS ";
                fprintf(fid,"%s\n",head_string);
                num_wan=0;
                for i=1:TBkitobj.basis_num
                    tempnum=maprule4(Projector_list(i));
                    num_wan = num_wan+tempnum;
                    fprintf(fid,"%d ",tempnum);
                end
                fprintf(fid,"\n");
                for i=1:TBkitobj.basis_num
                    if Projector_list(i) ~= -1
                        fprintf(fid,"%s ",'X');
                        % fprintf(fid,"%s ",TBkitobj.elementL(i)); % here we can
                        fprintf(fid,"%s ",maprule3(Projector_list(i)));
                        fprintf(fid,"\n");
                    end
                end
                fclose(fid);
            end
        end
