function gen_list = read_Magnetic_Sym_El(filename, crystal_system, options)
% read generators of a magnetic group from Mvasp2trace format file
arguments
    filename
    crystal_system {mustBeMember(crystal_system, [ ...
        "cubic", ...
        "tetragonal", ...
        "orthorhombic", ...
        "hexagonal", ...
        "trigonal", ...
        "triclinic", ...
        "monoclinic", ...
        ])}
    options.read_translation logical = false
end
if isa(filename,'char') || isa(filename,'string')
    fid = fopen(filename);
    line_1 = fgetl(fid);
    oper_num = str2double(line_1);
    FileMode = true;
elseif ~isscalar(filename)
    load("Table_of_Magnetic_Sym.mat");
    MSG_number1 = filename(1);
    MSG_number2 = filename(2);
    id = Table_of_Magnetic_Sym.mag_space_group == MSG_number1 & Table_of_Magnetic_Sym.number2 == MSG_number2;
    Table_of_Magnetic_Sym_check = Table_of_Magnetic_Sym(id,:);
    oper_num = Table_of_Magnetic_Sym_check.Number_of_Symmetry_operations;
    FileMode = false;
end

switch crystal_system
    case "cubic"
        R_matrix = eye(3);
    case "tetragonal"
        R_matrix = eye(3);
    case "orthorhombic"
        R_matrix = eye(3);
    case "hexagonal"
        R_matrix = [[1, -1/2, 0]; [0, sym(sqrt(3))/2, 0]; [0, 0, 1]];
    case "trigonal"
        R_matrix = eye(3);
        error("Not supported yet")
    case "triclinic"
        R_matrix = eye(3);
        error("Not supported yet")
    case "monoclinic"
        R_matrix = eye(3);
        error("Not supported yet")
end


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