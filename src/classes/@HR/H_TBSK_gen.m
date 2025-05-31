function H_hr = H_TBSK_gen(H_hr,options)
% H_TBSK_GEN Generate tight-binding Hamiltonian with Slater-Koster parameterization
%
%   H_hr = H_TBSK_GEN(H_hr,options) constructs a tight-binding Hamiltonian
%   using Slater-Koster parameterization with symbolic coefficients.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object with neighbor information
%       options - Structure with parameters:
%           level_cut: Maximum hopping level (default: -1)
%           onsite: Include onsite terms (logical, default: false)
%           per_dir: Periodic directions (default: [1,1,1])
%           chiral: Chiral term filtering (logical, default: false)
%           spin: Spin term filtering (logical, default: false)
%           SymHopping_Human: Custom symbolic hopping terms (optional)
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object with generated Hamiltonian
%
%   NOTES:
%       - Requires precomputed nn_store_smart
%       - Uses symbolic representation for coefficients
%       - Includes progress reporting
%
%   SEE ALSO:
%       HR, TBSK_Var_gen
%

arguments
    H_hr HR;
    options.level_cut {mustBeInteger} = -1;
    options.onsite logical = false;
    options.per_dir double = [1,1,1];
    options.chiral logical = false;
    options.spin logical = false;
    options.method  = 'nn_sparse';
    options.rough = false;
    options.vectorL = [];
    options.SymHopping_Human = [];
end
H_hr(1).num = false;
H_hr(1).coe = true;
if isempty(H_hr(1).nn_store)
    if isempty(H_hr(1).nn_store_smart)
        error('You have not run the function: nn/nn_smart');
    else
        if strcmp(options.method,'nn_sparse')
            error('You should use the key value: ''method'',''nn_smart'' as input.');
        end
    end
end
N_orbit = H_hr.WAN_NUM;
if options.level_cut >0
    select_nn_store = H_hr(1).nn_store(H_hr.nn_store(:,10)<=options.level_cut,:);
else
    select_nn_store = H_hr(1).nn_store;
end
if options.chiral
    i_L = select_nn_store(:,1);
    j_L = select_nn_store(:,2);
    element1_L = H_hr(1).elementL(i_L);
    element2_L =  H_hr(1).elementL(j_L);
    select_nn_store =select_nn_store(element1_L ~= element2_L,:);
end
if options.spin
    i_L = select_nn_store(:,1);
    j_L = select_nn_store(:,2);
    spin1_L = H_hr(1).quantumL(i_L,4);
    spin2_L = H_hr(1).quantumL(j_L,4);
    select_nn_store =select_nn_store(spin1_L == spin2_L,:);
end
Rvector_L = select_nn_store(:,6:8);
[Rvector_L_unique,sorted_label,cut_slice] = TBkit.cut_tools(Rvector_L);
select_nn_store = select_nn_store(sorted_label,:);
i_L = select_nn_store(:,1);
j_L = select_nn_store(:,2);
Rlength_L = select_nn_store(:,9);
l_L = select_nn_store(:,3)./Rlength_L;
m_L = select_nn_store(:,4)./Rlength_L;
n_L = select_nn_store(:,5)./Rlength_L;
nn_level_L = select_nn_store(:,10);
L_1_L = H_hr(1).quantumL(i_L,2);
L_2_L = H_hr(1).quantumL(j_L,2);
m_1_L = H_hr(1).quantumL(i_L,3);
m_2_L = H_hr(1).quantumL(j_L,3);
NRPTS_ = size(cut_slice,1);
fprintf('Generating Symbolic Hopping term ...')
if ~H_hr(1).overlap
    SymHopping_L = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L);
else
    SymHopping_L = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap',false);
    try
        SymHopping_L_str = string(SymHopping_L);
        SymHopping_L_str = strrep(SymHopping_L_str,"V","S");
        SymOverlap_L = str2sym(SymHopping_L_str);
    catch
        SymOverlap_L = HR.TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap',true);
        %error('')
    end
end
if ~isempty(options.SymHopping_Human)
    unSHL = unique(SymHopping_L);
    SymHopping_L = subs(SymHopping_L, unSHL, options.SymHopping_Human');
end
fprintf('< ok >\n')
if strcmp(H_hr.Type,'mat')
    H_hr(1) = H_hr(1).add_empty_one(Rvector_L_unique,'OverlapTest',false);   
end
pb = CmdLineProgressBar('Setting NRPT ');
if H_hr(1).overlap
    H_hr(2) = H_hr(1);
end
for i = 1:NRPTS_
    if ~H_hr(1).overlap
        H_hr = H_hr.set_hop(...
            SymHopping_L(cut_slice(i,1):cut_slice(i,2)),...
            i_L(cut_slice(i,1):cut_slice(i,2)),...
            j_L(cut_slice(i,1):cut_slice(i,2)),...
            Rvector_L_unique(i,:),'symadd');
    else
        H_hr(1) = H_hr(1).set_hop(...
            SymHopping_L(cut_slice(i,1):cut_slice(i,2)),...
            i_L(cut_slice(i,1):cut_slice(i,2)),...
            j_L(cut_slice(i,1):cut_slice(i,2)),...
            Rvector_L_unique(i,:),'sym'); 
        H_hr(2) = H_hr(2).set_hop(...
            SymOverlap_L(cut_slice(i,1):cut_slice(i,2)),...
            i_L(cut_slice(i,1):cut_slice(i,2)),...
            j_L(cut_slice(i,1):cut_slice(i,2)),...
            Rvector_L_unique(i,:),'sym');
    end
    pb.print(i,NRPTS_,' ...');
end
if options.onsite
    for i=1:N_orbit
        onsite_sym_name = "E__"+string(H_hr(1).elementL(i,1))+"_"...
            +string(H_hr(1).quantumL(i,2));...
            H_hr(1) = H_hr(1).set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'symadd');
    end
end
if H_hr(1).overlap
    for i=1:N_orbit
        H_hr(2) = H_hr(2).set_hop_single(1,i,i,[0 0 0],'symadd');
    end
end
pb.delete();
end
