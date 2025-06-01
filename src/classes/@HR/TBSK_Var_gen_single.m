function varargout = TBSK_Var_gen_single(L_1,L_2,m_1,m_2,nn_level,l,m,n,options)
%TBSK_VAR_GEN_SINGLE Generate Slater-Koster variables for single orbital pair
%
%   [COEFF, V_STR, S_STR] = TBSK_VAR_GEN_SINGLE(L_1, L_2, M_1, M_2, NN_LEVEL, L, M, N)
%   Generates Slater-Koster coefficients and variable names for a single orbital pair
%
%   Inputs:
%       L_1, L_2   - Angular momentum quantum numbers (0-3)
%       m_1, m_2   - Magnetic quantum numbers (-L to L)
%       nn_level   - Nearest neighbor level (default = -1)
%       l, m, n    - Direction cosines (default = 0)
%       options.overlap - Include overlap terms (default = false)
%
%   Outputs:
%       Coeff      - Slater-Koster coefficients [3×1]
%       V_str      - Hopping variable names {'VssS', 'VssP', 'VssD'}
%       S_str      - Overlap variable names (when options.overlap=true)
%
%   Note:
%       Automatically handles L_1 > L_2 case with proper sign changes
%
%   See also TBSK_VAR_GEN, HR.TBSK_Coeff_gen
arguments
    L_1 double{mustBeInteger};
    L_2 double{mustBeInteger};
    m_1 double{mustBeInteger};
    m_2 double{mustBeInteger};
    nn_level double{mustBeInteger} = -1;
    l  = 0;
    m  = 0;
    n  = 0;
    options.overlap = false;
    options.sym_mode = false;
end
if options.sym_mode || isa(l,'sym') || isa(m,'sym') ||isa(n,'sym')
    if ~isa(l,'sym')
        syms l real;
    end
    if ~isa(m,'sym')
        syms m real;
    end
    if ~isa(n,'sym')
        syms n real;
    end
end
optionscell = namedargs2cell(options);
if L_1 > L_2
    if nargout == 2
        [varargout{1},varargout{2} ] = HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,optionscell{:});
    else
        varargout{1} = HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,nn_level,-l,-m,-n,optionscell{:});
    end
    return;
end
if options.sym_mode
    Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n,'sym_mode',true);
else
    Coeff = HR.TBSK_Coeff_gen(L_1,L_2,m_1,m_2,l,m,n);
end
if options.overlap
    Char0 = 'S';
else
    Char0 = 'V';
end
switch L_1
    case 0
        Char1 = 's';
    case 1
        Char1 = 'p';
    case 2
        Char1 = 'd';
    case 3
        Char1 = 'f';
end
switch L_2
    case 0
        Char2 = 's';
    case 1
        Char2 = 'p';
    case 2
        Char2 = 'd';
    case 3
        Char2 = 'f';
end

E_hop{1} = [Char0,Char1,Char2,'S'] ;
E_hop{2} = [Char0,Char1,Char2,'P'] ;
E_hop{3} = [Char0,Char1,Char2,'D'] ;
if nn_level > 0
    E_hop{1} = [E_hop{1},'_',num2str(nn_level)];
    E_hop{2} = [E_hop{2},'_',num2str(nn_level)];
    E_hop{3} = [E_hop{3},'_',num2str(nn_level)];
end

if nargout == 2
    varargout{1} = Coeff;
    varargout{2} = E_hop;
else
    varargout{1} = Coeff;
end

end
