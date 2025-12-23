function gen_list = read_Magnetic_Sym_El(filename, MSG_number1, options)
% read generators of a magnetic group from Mvasp2trace format file
arguments
    filename
    MSG_number1 = 0
    options.read_translation logical = false
end
if isa(filename,'char') || isa(filename,'string')
    fid = fopen(filename);
    line_1 = fgetl(fid);
    oper_num = str2double(line_1);
    FileMode = true;
    if MSG_number1 == 0
        error("You must specify the space group number in order to judge the type of Bravais lattices")
    end
elseif ~isscalar(filename)
    load("Table_of_Magnetic_Sym.mat");
    MSG_number1 = filename(1);
    MSG_number2 = filename(2);
    id = Table_of_Magnetic_Sym.mag_space_group == MSG_number1 & Table_of_Magnetic_Sym.number2 == MSG_number2;
    Table_of_Magnetic_Sym_check = Table_of_Magnetic_Sym(id,:);
    oper_num = Table_of_Magnetic_Sym_check.Number_of_Symmetry_operations;
    FileMode = false;
end

%% 14 Bravais lattices, defined by The Mathematical Theory of Symmetry in Solids Representation Theory for Point Groups and Space Groups (Christopher Bradley, Arthur Cracknell)
if     ismember(MSG_number1, 1:2) % triclinic
    error("Not supported yet")
elseif ismember(MSG_number1, 3:15) % monoclinic
    error("Not supported yet")
elseif ismember(MSG_number1, [16:19,25:34,47:62]) % orthorhombic-P
    R_matrix = [[0, 1, 0]; [-1, 0, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [20:21,35:41,63:68]) % orthorhombic-C
    R_matrix = [[1/2, 1/2, 0]; [-1/2, 1/2, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [23:24,44:46,71:74]) % orthorhombic-I
    R_matrix = [[1/2, -1/2, 1/2]; [1/2, -1/2, -1/2]; [1/2, 1/2, -1/2]];
elseif ismember(MSG_number1, [22,42:43,69:70]) % orthorhombic-F
    R_matrix = [[1/2, 0, 1/2]; [0, -1/2, -1/2]; [1/2, 1/2, 0]];
elseif ismember(MSG_number1, [75:78,81,83:86,89:96,99:106,111:118,123:138]) % Tetragonal-P
    R_matrix = eye(3);
elseif ismember(MSG_number1, [79:80,82,87:88,97:98,107:110,119:122,139:142]) % Tetragonal-I  
    R_matrix = [[-1/2, 1/2, 1/2]; [1/2, -1/2, 1/2]; [1/2, 1/2, -1/2]];
elseif ismember(MSG_number1, 143:167) % trigonal
    R_matrix = [[0, sym(sqrt(3))/2, -sym(sqrt(3))/2]; [-1, 1/2, 1/2]; [1, 1, 1]];
elseif ismember(MSG_number1, 168:194) % hexagonal
    R_matrix = [[0, sym(sqrt(3))/2, 0]; [-1, 1/2, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [195,198,200,201,205,207,208,212,213,215,218,221:224]) % Cubic-P
    R_matrix = eye(3);
elseif ismember(MSG_number1, [196,202,203,209,210,216,219,225:228]) % Cubic-F
    R_matrix = [[0, 1/2, 1/2]; [1/2, 0, 1/2]; [1/2, 1/2, 0]];
elseif ismember(MSG_number1, [197,199,204,206,211,214,217,220,229,230]) % Cubic-I  
    R_matrix = [[-1/2, 1/2, 1/2]; [1/2, -1/2, 1/2]; [1/2, 1/2, -1/2]];
end

%%
gen_list = repmat(Oper(), 1, oper_num);
for i = 1:oper_num
    if FileMode
        line_i = fgetl(fid);
        oper_i_list_form = str2num(line_i); % must use str2num here
        oper_i_isUnitary = oper_i_list_form(end)==1;
    else
        oper_i_list_form = Table_of_Magnetic_Sym_check.symmetry_operation_R{1}(i,:);
        oper_i_isUnitary = Table_of_Magnetic_Sym_check.unitary_antiunitary{1}(i)==1;
    end
    oper_i_R = reshape(oper_i_list_form(1:9),[3,3])';
    oper_i_R = R_matrix * oper_i_R * R_matrix^(-1);

    if options.read_translation
        if FileMode
            oper_i_t = oper_i_list_form(10:12);

        else
            oper_i_t = Table_of_Magnetic_Sym_check.symmetry_operation_t{1}(i,:);
        end
    else
        oper_i_t = [0 0 0];
    end
    oper_i = Oper(double(oper_i_R), NaN, oper_i_t, "conjugate",~oper_i_isUnitary);
    gen_list(i) = oper_i;
end

if FileMode
    fclose(fid);
else
end

end