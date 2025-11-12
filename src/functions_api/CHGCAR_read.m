function [CHGCAR, MAG_DENSITY] = CHGCAR_read(filename)
% see https://www.vasp.at/wiki/index.php/CHGCAR for more detailed
% description of the CHGCAR file
arguments
    filename = "CHGCAR"
end
DENSITY_SET = CHGCAR_raw_read(filename);
CHGCAR = DENSITY_SET{1};
switch length(DENSITY_SET)
    case 1 % no-magnetic
        MAG_DENSITY = 0;
    case 2 % spin-polarized
        MAG_DENSITY = DENSITY_SET{2};
    case 4 % noncollinear
        MAG_DENSITY = cat(4, DENSITY_SET{2}, DENSITY_SET{3}, DENSITY_SET{4});
end
end
