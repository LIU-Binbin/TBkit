function [POT1, POT2] = LOCPOT_read(filename)
% see https://www.vasp.at/wiki/LOCPOT for more detailed
% description of the LOCPOT file
arguments
    filename = "LOCPOT"
end
DENSITY_SET = CHGCAR_raw_read(filename);
POT1 = DENSITY_SET{1};
switch length(DENSITY_SET)
    case 1 % no-magnetic: scalar potential
        POT2 = 0;
    case 2 % spin-polarized: spin up and spin down potentials
        POT2 = DENSITY_SET{2};
    case 4 % noncollinear: one scalar potential and three B-field-like potentials
        POT2 = cat(4, DENSITY_SET{2}, DENSITY_SET{3}, DENSITY_SET{4});
end
end
