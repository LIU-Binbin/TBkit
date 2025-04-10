function [coeff_trig, symvar_list_trig, H_htrig] = split_sym_eq(H_htrig, symbolic_polynomial)
%SPLIT_SYM_EQ  Splits a symbolic polynomial into coefficients and variable lists for symbolic Hamiltonians.
%
%   [coeff_trig, symvar_list_trig, H_htrig] = SPLIT_SYM_EQ(H_htrig, symbolic_polynomial)
%   splits the given symbolic polynomial into its coefficients and associated 
%   symbolic variables (momentum components) and updates the Htrig object (H_htrig).
%   This function processes the Hamiltonian terms, rewriting them in various forms 
%   (e.g., exponential, sine/cosine) and extracting the corresponding coefficients.
%
%   Inputs:
%       H_htrig            - An instance of the Htrig class, containing symbolic
%                             Hamiltonian trigonometric terms (HsymL_trig_bk).
%       symbolic_polynomial - The symbolic polynomial to be split. This represents
%                             the Hamiltonian or interaction terms involving momentum.
%
%   Outputs:
%       coeff_trig          - The coefficients of the terms in the symbolic polynomial.
%                             These coefficients are the numerical or symbolic parts of
%                             the Hamiltonian terms that correspond to specific momentum components.
%       symvar_list_trig    - The list of symbolic variables (e.g., k_x, k_y, k_z) involved 
%                             in the polynomial, representing the momentum components in the Hamiltonian.
%       H_htrig             - Updated Htrig object, potentially with new symbolic variables or 
%                             coefficients added based on the symbolic polynomial.
%
%   Behavior:
%       - Depending on the type of the Hamiltonian (`'exp'`, `'mat'`, `'list'`, `'sincos'`), 
%         the polynomial is rewritten and simplified into the corresponding form.
%       - The momentum components (k_x, k_y, k_z, etc.) are standardized across different
%         conventions used in the symbolic polynomial.
%       - The coefficients and the symbolic variables are extracted from the rewritten polynomial.
%       - If needed, additional basis functions may be added to the Htrig object.
%
%   Example:
%       syms kx ky kz real
%       H_htrig = Htrig();              % Create symbolic Htrig object
%       symbolic_polynomial = kx^2 + ky^2 + kz^2; % Example Hamiltonian term
%       [coeff_trig, symvar_list_trig, H_htrig] = split_sym_eq(H_htrig, symbolic_polynomial);
%
%   See also: Htrig, coeffs, simplify, rewrite, combine

    switch H_htrig.Type
        case {'exp'}
            k_symbol = combine(rewrite(symbolic_polynomial,'exp'));
        case {'mat','list'}
            k_symbol = expand(simplify(rewrite(symbolic_polynomial,'exp'),'IgnoreAnalyticConstraints',true));
        case 'sincos'
            k_symbol = expand(simplify(rewrite(symbolic_polynomial,'sincos'),'IgnoreAnalyticConstraints',true));
        otherwise
            k_symbol = expand(rewrite(simplify(symbolic_polynomial),'sincos'));
    end
    k_symbol_str = string(k_symbol);
    symvar_list = symvar(k_symbol);
    for i = 1:length(symvar_list)
        str_tmp = string(symvar_list(i));
        switch str_tmp
            case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'}
                k_symbol_str = strrep(k_symbol_str,str_tmp,'k_x');
            case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'}
                k_symbol_str = strrep(k_symbol_str,str_tmp,'k_y');
            case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'}
                k_symbol_str = strrep(k_symbol_str,str_tmp,'k_z');
            case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
                k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
        end
    end
    switch H_htrig.Type
        case {'exp','mat','list'}
            k_symbol = combine(str2sym(k_symbol_str),'exp');
        otherwise
            k_symbol = expand(combine(str2sym(k_symbol_str),'sincos'));
    end
    k_symbol_children1 = children(k_symbol);
    if isequal(simplify(fold(@mtimes,[k_symbol_children1{:}])),k_symbol)
        k_symbol_children{1} = k_symbol;
    elseif isequal(simplify(fold(@plus,[k_symbol_children1{:}])),k_symbol)
        k_symbol_children = k_symbol_children1;
    else
        k_symbol_children{1} = k_symbol;
    end
    switch H_htrig.Type
        case {'mat','list'}
            nSon = length(k_symbol_children);
            symvar_list_trig = zeros(nSon,3,'sym');
            coeff_trig = zeros(nSon,1,'sym');
            if k_symbol == sym(0)
            else
                for k = 1:nSon
                    [TmpCoe,TmpVar] = coeffs(k_symbol_children{k});
                    ActualVar = simplify(log(TmpVar),'IgnoreAnalyticConstraints',true)/1i;
                    [ActualCoe,ktype] = coeffs(ActualVar,H_htrig.seedsvar);
                    if ismember(sym(1),ktype)
                        n1 = find(sym(1)==ktype);
                        TmpCoe = simplify(exp(ActualCoe(n1)*(1i)),'IgnoreAnalyticConstraints',true)*TmpCoe;
                        ActualCoe(n1) = [];
                        ktype(n1) = [];
                    end
                    HsymL_coeLtmp = zeros(1,3,'sym');
                    for i = 1:3
                        ik = find(H_htrig.seedsvar(i) == ktype);
                        if ~isempty(ik)
                            HsymL_coeLtmp(i) = ActualCoe(ik);
                        end
                    end
                    symvar_list_trig(k,:) = HsymL_coeLtmp;
                    coeff_trig(k) = TmpCoe;
                end
            end
            if strcmp(H_htrig.Type,'mat')
                H_htrig = H_htrig.add_empty_one(symvar_list_trig);
            end
            return;
        otherwise
            for k = 1:length(k_symbol_children)
                try
                    [coeff_trig,~] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
                catch
                    coeff_trig = k_symbol_children{k};
                end
                for i = 1:numel(coeff_trig)
                    tmp_label = contains(string(coeff_trig),H_htrig.seeds);
                    if sum(tmp_label)
                        [~,coeff_trig_list] = coeffs(coeff_trig(i));
                        for j = 1:numel(coeff_trig_list)
                            tmp_label2 = contains(string(coeff_trig_list(j)),H_htrig.seeds);
                            if sum(tmp_label2)
                                H_htrig = H_htrig.find_HsymL_trig_bk(coeff_trig_list(j));
                            end
                        end
                    end
                end
            end
    end
    coeff_trig = sym([]);
    symvar_list_trig= sym([]);
    for k = 1:length(k_symbol_children)
        [coeff_trig_tmp,symvar_list_trig_tmp] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
        coeff_trig = [coeff_trig,coeff_trig_tmp];
        symvar_list_trig = [symvar_list_trig,symvar_list_trig_tmp];
    end
end

