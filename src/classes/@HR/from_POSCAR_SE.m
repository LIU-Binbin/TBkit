function H_hr = from_POSCAR_SE(POSCAR_file,options)
% FROM_POSCAR_SE Create HR object from POSCAR file
%
%   H_hr = FROM_POSCAR_SE(POSCAR_file,options) generates an HR object
%   from a POSCAR file with various configuration options.
%
%   INPUT ARGUMENTS:
%       POSCAR_file - Path to POSCAR file (default: 'POSCAR')
%       options - Structure with parameters:
%           WAN_NUM: Number of Wannier orbitals (-1 for auto)
%           r_max: Maximum hopping distance
%           level_cut: Energy level cutoff
%           per_dir: Periodic directions [1,1,1]
%           chiral: Include chiral terms (logical)
%           spin: Include spin terms (logical)
%           deltarule: Delta rule type (0-2)
%           alpharule: Alpha rule type (0-2)
%           onsite: Include onsite terms (logical)
%           Type: Output type ('mat','list','sparse')
%           overlap: Include overlap matrix (logical)
%           symbolic: Use symbolic representation (logical)
%           E0: Onsite energy shift
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object generated from POSCAR
%
%   NOTES:
%       - Uses TBkit.POSCAR_read for file parsing
%       - Supports different Hamiltonian generation rules
%
%   SEE ALSO:
%       HR, TBkit.POSCAR_read
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
    POSCAR_file char= 'POSCAR';
    options.WAN_NUM double{mustBeInteger}= -1;
    options.r_max = 3.3;
    options.level_cut = -1;
    options.per_dir double = [1,1,1];
    options.chiral logical = false;
    options.spin logical = false;
    options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
    options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
    options.onsite logical = true;
    options.Accuracy = 1e-3;
    options.search_range = [2 2 2];
    options.Type char{mustBeMember(options.Type,{'mat','list','sparse'})}= 'mat';
    options.overlap = false;
    options.symbolic = false;
    options.simplify = false;
    options.Rd  = -1;
    options.para struct=struct();
    options.silence = true;
    options.E0 = 0;
end
r_max_search =options.r_max ;
level_cut  = options.level_cut;
per_dir = options.per_dir;
onsite = options.onsite;
Accuracy = options.Accuracy ;
search_range = options.search_range ;
Type = options.Type;
[~,sites,~,~,~]=TBkit.POSCAR_read(POSCAR_file,'tbsk');
if options.WAN_NUM == -1
    H_hr = HR(length(sites),'overlap',options.overlap,'Type',Type);
else
    H_hr = HR(options.WAN_NUM,'overlap',options.overlap,'Type',Type);
    %H_hr(1,2) = H_hr;
end
H_hr.Basis_num = H_hr.WAN_NUM;
if options.symbolic
    H_hr =  POSCAR_file > H_hr;
else
    H_hr = H_hr < POSCAR_file;
end
switch Type
    case 'mat'
        H_hr.Type = 'mat';
        H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
        H_hr = H_hr.H_TBSK_gen(...
            'chiral', options.chiral,...
            'spin',options.spin ,...
            'onsite',logical(onsite),...
            'per_dir',per_dir ,...
            'level_cut',level_cut ...
            );
        if isequal(sym(options.E0),(0))

        else
            H_hr(1) = H_hr(1).set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'symadd');
        end
    case 'list'
        H_hr.Type = 'list';
        H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
        H_hr = H_hr.H_TBSK_gen(...
            'chiral', options.chiral,...
            'spin',options.spin ,...
            'onsite',logical(onsite),...
            'per_dir',per_dir ,...
            'level_cut',level_cut ...
            );
        if isequal(sym(options.E0),(0))
        else
            H_hr(1) = H_hr(1).set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'symadd');
        end
    case 'sparse'
        H_hr = H_hr.nn(search_range,Accuracy,r_max_search);
        H_hr = H_hr.H_TBSK_gen_sparse(...
            'chiral', options.chiral,...
            'spin',options.spin ,...
            'onsite',logical(onsite),...
            'per_dir',per_dir ,...
            'level_cut',level_cut,...
            'Rd',options.Rd,...
            'para',options.para,...
            'deltarule',options.deltarule,...
            'alpharule',options.alpharule ...
            );
        if isequal(sym(options.E0),(0))
        else
            H_hr(1) = H_hr(1).set_hop_mat(eye(H_hr.WAN_NUM)*options.E0,[0,0,0],'add');
        end
end
if ~strcmp(Type,'sparse')
    if options.deltarule
        fprintf('applying delta rule ...\n');
        H_hr = H_hr.deltarule(level_cut,options.deltarule,'Rd',options.Rd);
    end
    if options.alpharule
        fprintf('applying alpha rule ...\n');
        H_hr = H_hr.alpharule(level_cut,options.alpharule,'Rd',options.Rd,'silence',options.silence);
    end
    if options.simplify
        fprintf('simplify the Hamitoian...\n');
        H_hr =H_hr.simplify();
    end
    disp(H_hr.symvar_list);
end
end
