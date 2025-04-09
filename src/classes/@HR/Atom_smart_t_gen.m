function Atom_smart_t = Atom_smart_t_gen(site1,site2)
%ATOM_SMART_T_GEN Generate smart atom transition descriptor
%   Creates structured descriptor for atomic transitions between two sites
%   Contains fractional coordinates, orbital information, and transition metadata
%
%   Inputs:
%       site1   - Source site structure with fields:
%                 rc1, rc2, rc3  : Fractional coordinates (3 components)
%                 seq            : Sequence number
%                 orb            : Orbital name
%                 orb_sym        : Orbital symmetry
%                 name           : Site name
%       site2   - Target site structure (same fields as site1)
%
%   Output:
%       Atom_smart_t - Transition descriptor structure with fields:
%           R_fractional_from  : Source site coordinates [rc1, rc2, rc3]
%           R_fractional_to    : Target site coordinates [rc1, rc2, rc3]
%           R_fractional_diff  : Coordinate difference (to - from)
%           seq_from           : Source sequence ID
%           seq_to             : Target sequence ID
%           l_name_from        : Source orbital name
%           l_name_to          : Target orbital name
%           orb_sym_from       : Source orbital symmetry
%           orb_sym_to         : Target orbital symmetry
%           handyname          : Formatted transition name "from -> to"

% Extract fractional coordinates from both sites
Rc1 = [site1.rc1, site1.rc2, site1.rc3];
Rc2 = [site2.rc1, site2.rc2, site2.rc3];

% Populate transition descriptor structure
Atom_smart_t = struct(...
    'R_fractional_from', Rc1,...
    'R_fractional_to',   Rc2,...
    'R_fractional_diff', Rc2 - Rc1,...  % Direct to-from calculation
    'seq_from',     site1.seq,...
    'seq_to',       site2.seq,...
    'l_name_from',  site1.orb,...
    'l_name_to',    site2.orb,...
    'orb_sym_from', site1.orb_sym,...
    'orb_sym_to',   site2.orb_sym,...
    'handyname',    strcat(site1.name,' -> ',site2.name)...
    );
end