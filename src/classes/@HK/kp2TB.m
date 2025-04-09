function varargout = kp2TB(H_hk,kpoints_f,groups,options)
%KP2TB Convert k·p Hamiltonian to tight-binding model
%
% Syntax:
%   HR = kp2TB(H_hk)
%   [HR,Htrig] = kp2TB(H_hk,kpoints_f,groups,options)
%
% Inputs:
%   H_hk - HK k·p Hamiltonian
%   kpoints_f - Folding k-points (default=[0,0,0])
%   groups - Symmetry groups for constraints
%   options - Optional parameters:
%       Accuracy - Numerical tolerance (default=1e-4)
%       search_range - Neighbor search range (default=[1 1 1])
%       R_max - Maximum hopping distance (default=-1: auto)
%       level_cut - Cutoff level for hoppings (default=1)
%       mini - Minimal model generation (default=false)
%
% Outputs:
%   HR - Resulting tight-binding HR object
%   Htrig - Trigonometric Hamiltonian (optional)
%
% Description:
%   Converts continuous k·p model to discrete tight-binding model by:
%   1. Mapping k-space terms to real-space hoppings
%   2. Applying symmetry constraints if provided
%   3. Solving equation systems for parameter matching
%
% Supported Modes:
%   - 'Gamma': Γ-point folding
%   - 'kpoints': General k-point folding
%   - 'TBsubs': Tight-binding with substitutions
%
% Example:
%   H_tb = kp2TB(H_kp,[0 0 0],symm_group,'Accuracy',1e-5)
arguments
    H_hk HK;
    kpoints_f double= [0,0,0];
    groups = [];
    options.Accuracy = 1e-4;
    options.search_range = [1 1 1];
    options.R_max = -1;
    options.level_cut = 1;
    options.mini = false;
end
if isempty(groups)
    if isequal(kpoints_f,[0,0,0])
        mode = 'Gamma';
    else
        mode = 'kpoints';
    end
else
    mode = 'TBsubs';
end
R_struct = H_hk.Rm2abc(H_hk.Rm);
if ~strcmp(mode,'TBsubs')
    if abs(abs(cosd(R_struct.gamma))) < 1e-4
        lattice_mode = 'T';
    elseif abs((cosd(R_struct.gamma)) - -1/2) < 1e-4
        lattice_mode = 'H';
    elseif abs((cosd(R_struct.gamma)) - 1/2) < 1e-4
        lattice_mode = 'H2';
    else
        error('unsupport yet! only tetragonal of hexagonal\n');
    end
    if strcmp(mode , 'Gamma')
    else
        warning('This function is not well-written, only support Term mode now');
        kpoints = kpoints_f * H_hk.Gk;
        H_hk = HK(H_hk.Basis_num,H_hk.Degree,H_hk.Term_to_save.subsk(kpoints));
    end
end
H_hr = HR(H_hk.Basis_num,'sym',true);
H_hr = H_hr.TBkitCopy(H_hk);
H_hr.coe = true;
H_hr.num = false;
if strcmp(mode,'TBsubs')
    if options.R_max == -1
        R_a = norm(H_hr.Rm(1,:));
        R_b = norm(H_hr.Rm(2,:));
        R_c = norm(H_hr.Rm(3,:));
        R_max = max(options.search_range)*max([R_a,R_b,R_c])+options.Accuracy;
    else
        R_max = options.R_max;
    end
    Accuracy = options.Accuracy;
    NAccuracy = round(-log(Accuracy)/log(10));
    H_hr= H_hr.nn(options.search_range,Accuracy,R_max);
    H_hr = H_hr.init('level_cut',options.level_cut,"onsite",1);
    basis_num = H_hk.Basis_num;
    if options.mini
        nonzeromat = zeros(H_hk.Basis_num);
        for i = 1:length(groups)
            H_hk_bk = H_hk.applyU(groups(i)) ;
            nonzeromat =nonzeromat+sum(logical((H_hk_bk.HcoeL~=sym(0))),3);
        end
        if strcmp(H_hr.type,'list')
            H_hr = H_hr.rewind();
        end
        zeromat = (nonzeromat == 0);
        for i = 1:basis_num
            for j = 1:basis_num
                if zeromat(i,j)
                    H_hr.HcoeL(i,j,:) = sym(0);
                end
            end
        end
        H_hr = H_hr.rewrite();
    end
    H_hr = H_hr.applyOper(groups,'generator',true);
    [H_hr,Sublist,Unique_term] = H_hr.unique();
    fprintf('******** Simplify Hamiontonian ********\n');
    fprintf('----------   SymVarNum: %d   ----------\n',length(H_hr.symvar_list));
    H_htrig = H_hr.HR2Htrig();
    H_hr_hk = H_htrig.Htrig2HK(kpoints_f,'Order',H_hk.Degree);
    label_list = find(H_hr_hk.HcoeL~=sym(0));
    HcoeL_tmp_1 = (H_hr_hk.HcoeL(label_list));
    HcoeL_tmp_2 = (H_hk.HcoeL(label_list));
    Symvar_for_sub = symvar(HcoeL_tmp_1);
    Equationlist_r = real(HcoeL_tmp_1)==real(HcoeL_tmp_2);
    Equationlist_i = imag(HcoeL_tmp_1)==imag(HcoeL_tmp_2);
    Equationlist = vpa([Equationlist_r;Equationlist_i],NAccuracy);
    Equationlist = HK.isolateAll(Equationlist,Symvar_for_sub);
    Equationlist = unique(simplify(Equationlist));
    Solved = solve(Equationlist,Symvar_for_sub);
    Solved_value = sym(ones(1,length(Symvar_for_sub)));
    count = 0;
    for symvartmp = Symvar_for_sub
        count = count+1;
        if isempty(Solved.(char(symvartmp)))
            Solved_value(count) = sym(0);
        else
            Solved_value(count) = Solved.(char(symvartmp));
        end
    end
    HcoeLtmp = subs(H_hr.HcoeL,Symvar_for_sub,Solved_value);
    Allzero =length(symvar(HcoeLtmp));
    if Allzero == 0
        warning('No subs result, check by yourselt')
        Problem_hr.HR = H_hr;
        Problem_hr.Equationlist = Equationlist;
        Problem_hr.Sublist = Sublist;
        Problem_hr.Unique_term = Unique_term;
        varargout{1} = Problem_hr;
        if nargout == 2
            varargout{2} = H_htrig;
        end
        return;
    else
        H_hr.HcoeL = HcoeLtmp;
    end
    H_hr = H_hr.rewrite();
    H_hr.HcoeL = TBkit.cleanVar(H_hr.HcoeL,Accuracy);
    varargout{1} = H_hr;
    if nargout == 2
        H_htrig.HcoeL  = subs(H_htrig.HcoeL,Symvar_for_sub,Solved_value);
        varargout{2} = H_htrig;
    end
else
    for i = 1:H_hk.Kinds
        if  any(H_hk.HcoeL(:,:,i)~=sym(0),'all')
            [vector_list,Coeffs_list] = HK.HstrL_classify(H_hk.HstrL(i),R_struct,lattice_mode);
            for j =1:length(Coeffs_list)
                H_hr = H_hr.set_hop_mat(H_hk.HcoeL(:,:,i)*Coeffs_list(j),vector_list(j,:),'symadd');
            end
        end
    end
    varargout{1} = H_hr;
end
end
