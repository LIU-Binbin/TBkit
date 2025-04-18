classdef Htrig < TBkit & matlab.mixin.CustomDisplay
    properties
        HcoeL;
        HnumL;
        HsymL_trig = sym(1);
    end
    properties
        HsymL_coeL ;
        HsymL_numL ;
        Htrig_num   ;
        Nslab = [0,0,0];
        Htrig_latex ;
        HsymL_trig_bk ;
        seeds = [];
        seedsvar = sym([]);
        Sigmas        = []     ;
        rm_list                ;
        Trig_list =    Trig()  ;
        Type {mustBeMember(Type,{'mat','list','sincos','exp','slab','sparse'})}= 'sincos';
        Hmat_pre;
        symvarL;
    end
    properties(Dependent = true)
        Htrig_sym   ;
        HsymL;
        Kinds;
    end
    properties
        num = false;
        coe = true;
        Duality_vector_dist;
    end
    methods (Access = protected)
        propgrp = getPropertyGroups(~)
    end
    %% constuction
    methods
        function H_htrig = Htrig(BASIS_NUM,Trig_list,options,propArgs)
            arguments
                BASIS_NUM = 4;
                Trig_list =[];
                options.Type = 'sincos';
                propArgs.?TBkit;
            end
            %
            propArgsCell = namedargs2cell(propArgs);
            H_htrig = H_htrig@TBkit(propArgsCell{:});
            VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
            H_htrig.seedsvar = VarUsing;
            H_htrig.seeds = string(H_htrig.seedsvar);
            %syms k_x k_y k_z real;
            switch nargin
                case 1
                    if isnumeric(BASIS_NUM)
                        H_htrig.Basis_num = BASIS_NUM;
                        H_htrig.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                        H_htrig.HnumL = (zeros(BASIS_NUM,BASIS_NUM,1));
                        H_htrig.Trig_list = Trig()  ;
                        switch options.Type
                            case 'sincos'
                                H_htrig.HsymL_trig_bk = [sin(VarUsing) cos(VarUsing) ];
                            case 'exp'
                                H_htrig.HsymL_trig_bk = [exp(1i*VarUsing) exp(-1i*VarUsing)];
                            case 'mat'
                                H_htrig.HsymL_coeL = sym(zeros(1,H_htrig.Dim)); % k_x k_y k_z ...
                                H_htrig.HsymL_numL = zeros(1,H_htrig.Dim);
                            case 'list'
                                H_htrig.HsymL_coeL = [zeros(1,H_htrig.Dim,'sym'),ones(1,2,'sym')];
                                H_htrig.HcoeL = sym(0);
                                H_htrig.HsymL_numL = [zeros(1,H_htrig.Dim,'double'),ones(1,2,'double')];
                                H_htrig.HnumL = 0;
                        end
                        H_htrig.Type = options.Type;
                    else
                        if isa(BASIS_NUM,'sym')
                            Hsym = BASIS_NUM;
                            if size(Hsym,1) ~= size(Hsym,2)
                                error('square sym mat required!')
                            end
                            BASIS_NUM = size(Hsym,1);
                            H_htrig.Basis_num = BASIS_NUM;
                            H_htrig.HcoeL = sym(zeros(BASIS_NUM,BASIS_NUM,1));
                            H_htrig.HnumL = (zeros(BASIS_NUM,BASIS_NUM,1));
                            H_htrig.Trig_list = Trig()  ;
                            H_htrig = H_htrig + Hsym;
                        end
                    end
                case 2
                    opitonsCell = namedargs2cell(options);
                    H_htrig = Htrig(BASIS_NUM,opitonsCell{:},propArgsCell{:});
                    H_htrig = H_htrig + Trig_list;
            end
        end
    end
     %% get methods
    methods
        function symvarL = get.symvarL(H_htrig)
            symvarL = [H_htrig.symvar_list,symvar(H_htrig.HsymL)];
        end
        function Htrig_sym = get.Htrig_sym(H_htrig)
            Htrig_sym = sym(zeros(H_htrig.Basis_num,H_htrig.Basis_num));
            if H_htrig.num
                H_htrig.HcoeL = sym(H_htrig.HnumL);
                H_htrig.HsymL_coeL = sym(H_htrig.HsymL_numL);
                %H_htrig.HsymL_trig = sym(H_htrig.HsymL_trig);
            end
            switch H_htrig.Type
                case {'exp','sincos'}
                    try
                        for i =1:H_htrig.Kinds
                            Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*H_htrig.HsymL_trig(i);
                        end
                    catch
                        Htrig_sym =sym([]);
                    end
                case {'slab'}
                    try
                        for i =1:H_htrig.Kinds
                            Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*H_htrig.HsymL_trig(i);
                        end
                    catch
                        Htrig_sym =sym([]);
                    end
                    if strcmp(H_htrig.Type , 'slab' )
                        tmpHsymL_trig = H_htrig.HsymL_trig;
                        for i = 1:length(tmpHsymL_trig)
                            A = char(tmpHsymL_trig(i));
                            A(end)= A(end-3);
                            A(end-3)='N';
                            tmpHsymL_trig(i) = str2sym(A);
                        end
                        Htrig_sym = tril(Htrig_sym)+subs(triu(Htrig_sym,1),H_htrig.HsymL_trig,tmpHsymL_trig);
                    end
                case 'mat'
                    syms k_x k_y k_z;
                    HvarL = exp(1i*H_htrig.HsymL_coeL*[k_x;k_y;k_z]);
                    for i =1:H_htrig.Kinds
                        Htrig_sym = Htrig_sym + H_htrig.HcoeL(:,:,i)*HvarL(i);
                    end
                case 'list'
                    syms k_x k_y k_z;
                    Factorlist = H_htrig.HcoeL.*exp(1i*H_htrig.HsymL_coeL(:,1:3)*[k_x;k_y;k_z]);
                    [ij_unique,sumFactorlist] = HollowKnight.generalcontractrow(double(H_htrig.HsymL_coeL(:,4:5)),Factorlist);
                    Htrig_sym = zeros(H_htrig.Basis_num,'sym');
                    indL = sub2ind([H_htrig.Basis_num,H_htrig.Basis_num],ij_unique(:,1),ij_unique(:,2));
                    Htrig_sym(indL) = sumFactorlist;
            end

            
        end
        function Htrig_latex = get.Htrig_latex(H_htrig)
            Htrig_latex = latex(H_htrig.Htrig_sym);
        end
        function HsymL = get.HsymL(H_htrig)
            switch H_htrig.Type
                case {'exp','sincos'}
                    HsymL = H_htrig.HsymL_trig;
                case {'mat','list'}
                    if H_htrig.coe
                        HsymL = H_htrig.HsymL_coeL;
                    else
                        HsymL = H_htrig.HsymL_numL;
                    end
                otherwise
                    HsymL = H_htrig.HsymL_trig;
            end
        end
        function Kinds = get.Kinds(H_htrig)
            switch H_htrig.Type
                case {'exp','sincos'}
                    Kinds = length(H_htrig.HsymL_trig);
                case {'mat','list'}
                    Kinds = max(size(H_htrig.HsymL_numL,1),size(H_htrig.HsymL_coeL,1));
                otherwise
                    Kinds = length(H_htrig.HsymL_trig);
            end
        end
    end
    
    %%
    methods
        H_htrig = setup_rough(H_htrig,symbolic_polynomial,pauli_mat,silence)
        H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence)
        H_htrig = set_hop(H_htrig,SymHopping,SymVar,Indij)
        [coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,symbolic_polynomial)
        H_htrig = find_HsymL_trig_bk(H_htrig,coeff_trig)
        Kind = k_symbol2Kind(H_htrig,k_symbol)
        H_htrig = reseq(H_htrig,wan_list,kinds_list)
        H_htrig = add_empty_one(H_htrig,vector)
    end
    methods
        H_htrig2 = rewrite(H_htrig,mode)
        H_htrig = diff(H_htrig,dir,options)
        H_htrig = simplify(H_htrig,Accuracy,options)
    end
    methods
        [HoutL] = HCAR_gen(H_htrig,klist,options)
        varargout = diff_klist(H_htrig,dir,klist,options)
    end
    methods
        H_htrig = init(H_htrig,HsymL_trig_in,options)
        H_htrig = applyOper(H_htrig,SymOper,options)
        H_htrig = dualize(H_htrig)
        H_htrig = hermitize(H_htrig)
        H_htrig_bk = subsOper(H_htrig,SymOper)
        H_htrig = applyU(H_htrig,U,conjugate ,antisymmetry )
        [H_htrig,H_htrig2] = applyR(H_htrig,R)
        [H_htrig,H_htrig2] = applyRt(H_htrig,Rt)
        [H_htrig,Smat] = Smatgen(H_htrig,R,Accuracy)
        Factorlist_parity = factorlist_parity(H_htrig)
        H_htrig = nn(H_htrig,search_range,Accuracy,Rlength_cut)
    end
    methods
        C = plus(A,B)
        H_htrig = uminus(H_htrig)
        C = minus(A,B)
        C = mtimes(A,B)
        C = mtimes_inv(B,A)
        C = gt(B,A)
        C = lt(A,B)
        C = le(A,B)
        C = horzcat(A,B)
        H_htrig = ctranspose(H_htrig)
        H_htrig = conj(H_htrig)
        Htrig_sym = sym(H_htrig,options)
        Htrig_latex = latex(H_htrig)
    end
    methods
        H_htrig = translate(H_htrig,U)
        [H_htrig,EQL] = subsVar(H_htrig,varargin)
        H_htrig2 = subs(H_htrig,varargin)
        H_htrig = rotation(H_htrig,Rotation)
        H_htrig2 = discretize(H_htrig,Nslab,options)
        H_hr = Htrig2HR(H_htrig,options)
        H_hk = Htrig2HK(H_htrig,kpoints_f,options)
    end
    methods
        varargout = EIGENCAR_gen(H_htrig,options)
        [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_slab(H_htrig,options)
        [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_wire(H_htrig,options)
        H_htrig = Subsall(H_htrig,mode,AddtionVar)
        mat = HsymL_trig2mat(H_htrig,HsymL_trig)
        [num_label,coe_label,H_htrig] = NumOrCoe(H_htrig)
    end
    methods (Static)
        H_htrig = HR2Htrig(H_hr)
        Delta_Oper = Delta_Oper(Oper_list)
        sym_coe_list = coeff_extract(sym_term)
    end
end
