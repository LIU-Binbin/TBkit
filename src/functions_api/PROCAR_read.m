function [PROCAR_collection,EIGENCAR] = PROCAR_read(filename,opts)

        %     PROCAR_collection_data_type
        % ************************************************************************************************************************************************
        %     ----- ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot
        %           ion1                                .  
        %           ion2                                . 
        %             .       .     .      .       WEIGHTCAR      .       .      .     .
        %             .                                 .
        %             .                                 .
        % ************************************************************************************************************************************************
           


 
arguments
    filename = 'PROCAR';
    opts.SOC_flag = 0;
    opts.F_flag = 0;
end
    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('File not found');
    end
    fclose(fid);

    % Read and filter band data
    fileContent = fileread(filename);
    lines = splitlines(fileContent);

    % Remove comment lines and headers
    nonEmptyRows = ~cellfun(@(x) isempty(strtrim(x)), lines);
    validLines = ~contains(lines, {'#', 'ion','k-point','PROCAR','tot'});
    EIGENCARLines = contains(lines, {'energy'});

    dataLines = lines(validLines&nonEmptyRows);
    CommentLines = lines(~validLines);
    PROCARinformation = str2double(regexp(CommentLines{2}, '\d+', 'match'));
    %
    EIGENCARinformation = lines(EIGENCARLines);
    energy_pattern = '# energy\s+([-+]?\d+\.\d+)';
    energy_tokens = regexp(strjoin(EIGENCARinformation,' '), energy_pattern, 'tokens');
    energy_values = cellfun(@(x) str2double(x{1}), energy_tokens);
    %kpoints
    kpointsnum=PROCARinformation(1);
    %bands
    bandsnum=PROCARinformation(2);
    %ions
    ionsnum=PROCARinformation(3);
    %
    EIGENCAR  = reshape(energy_values, [bandsnum kpointsnum ]);
    %
    if opts.F_flag
        orbitalsnumPlus2=18;
    else
        orbitalsnumPlus2=11;
    end
    % Convert to numeric matrix
    if opts.SOC_flag
        ColNum  = ionsnum*kpointsnum*bandsnum*4;
    else
        ColNum  = ionsnum*kpointsnum*bandsnum;
    end
    all_data = strjoin(dataLines, ' ');
    data = sscanf(all_data, '%f', [orbitalsnumPlus2,ColNum ])';
    
    % WEIGHTCAR datatype
    % ************************************************************************************************************************************************
    %     ----- band     kpoint1     kpoint2     kpoint3     ............
    %           band1                                .
    %           band2                                .
    %             .       .     .      .      .   weight      .     .     .     .
    %             .                                  .
    %             .                                  .
    % ************************************************************************************************************************************************
    orbitalsnumPlus1 = orbitalsnumPlus2-1;
    if opts.SOC_flag
        PROCAR_collection = cell(ionsnum,orbitalsnumPlus1,4);
        for iorb = 1:orbitalsnumPlus1
            TempData = data(:,iorb+1);
            for ion = 1:ionsnum
                for i = 1:4
                    PROCAR_collection{ion,iorb,i} = reshape( ...
                        TempData((ion+i-1):(ionsnum*4):ColNum,:), ...
                        [bandsnum,kpointsnum]);
                end
            end
        end
    else
        PROCAR_collection = cell(ionsnum,orbitalsnumPlus1);
        for iorb = 1:orbitalsnumPlus1
            TempData = data(:,iorb+1);
            for ion = 1:ionsnum
                PROCAR_collection{ion,iorb} = reshape(TempData(ion:ionsnum:ColNum,:),[bandsnum,kpointsnum]);
            end
        end
    end

end
         
        
       
    





