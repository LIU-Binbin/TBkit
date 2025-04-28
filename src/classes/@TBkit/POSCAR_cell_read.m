function POSCAR = POSCAR_cell_read(filename,formatSpec)
            startRow = 2;
            delimiter = {'\t',' '};
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
            end
            rawStringColumns = string(raw(:,:));
            POSCAR=rawStringColumns;
            clearvars filename delimiter startRow formatSpec fileID dataArray ans rawStringColumns raw;
        end
