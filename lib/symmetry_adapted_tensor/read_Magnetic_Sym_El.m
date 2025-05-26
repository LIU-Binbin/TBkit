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

fid = fopen(filename);

line_1 = fgetl(fid);
oper_num = str2double(line_1);

hex2cart = [[1, -1/2, 0]; [0, sym(sqrt(3))/2, 0]; [0, 0, 1]];

gen_list = repmat(Oper(), 1, oper_num);
for i = 1:oper_num
    line_i = fgetl(fid);
    oper_i_list_form = str2num(line_i); % must use str2num here
    oper_i_R = reshape(oper_i_list_form(1:9),[3,3])';
    
    switch crystal_system
        case "cubic"
            true
        case "tetragonal"
            true
        case "orthorhombic"
            true
        case "hexagonal"
            oper_i_R = hex2cart * oper_i_R * hex2cart^(-1);
        case "trigonal"
            error("Not supported yet")
        case "triclinic"
            error("Not supported yet")
        case "monoclinic"
            error("Not supported yet")
    end

    if options.read_translation
        oper_i_t = oper_i_list_form(10:12);
    else
        oper_i_t = [0 0 0];
    end
    
    oper_i_isUnitary = oper_i_list_form(end)==1;
    
    oper_i = Oper(double(oper_i_R), NaN, oper_i_t, "conjugate",~oper_i_isUnitary);
    gen_list(i) = oper_i;
end
fclose(fid);

end