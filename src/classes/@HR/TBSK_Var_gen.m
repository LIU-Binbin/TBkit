function varargout = TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,options)
%TBSK_VAR_GEN Generate Tight-Binding Slater-Koster variables
%
%   [Coeff, V_str, S_str] = TBSK_VAR_GEN(L_1_L, L_2_L, m_1_L, m_2_L, nn_level_L, l_L, m_L, n_L)
%   Generates Slater-Koster coefficients and variable names for tight-binding calculations
%
%   Inputs:
%       L_1_L    - Angular momentum quantum numbers for first orbital (integer array)
%       L_2_L    - Angular momentum quantum numbers for second orbital (integer array)
%       m_1_L    - Magnetic quantum numbers for first orbital (integer array)
%       m_2_L    - Magnetic quantum numbers for second orbital (integer array)
%       nn_level_L - Nearest neighbor level (integer, default = -1)
%       l_L      - Direction cosine l component (default = 0)
%       m_L      - Direction cosine m component (default = 0)
%       n_L      - Direction cosine n component (default = 0)
%       options.overlap - Logical flag for overlap calculation (default = false)
%
%   Outputs:
%       Coeff    - Slater-Koster coefficients matrix
%       V_str    - Variable names for hopping terms
%   See also HR.TBSK_Var_gen_single, HR.TBSK_Coeff_gen
% <https://en.wikipedia.org/wiki/Tight_binding>
% L_1,L_2 : orbital angular momentum: s 0 p 1 d 2 f 3
% m_1,m_2 : magnetic quantum number:
%           - s 0
%           - p y:-1 z:0 x:1 % because py = i/sart(2)*(Y-1
%           +Y1), we label this i as a minus number
%           - d xy:-2 yz:-1 z2:0 xz:1 x2-y2:2
%           - f
%           y(3x2-y2):-3 xyz:-2  yz2:-1
%           z3:0
%           x(x2-3y2):3  z(x2-y2):2  xz2:1

% Coeff = zeros(nhopping,3);
% when generating a batch of Hoppings, we must code vectorized
arguments
    L_1_L double{mustBeInteger};
    L_2_L double{mustBeInteger};
    m_1_L double{mustBeInteger};
    m_2_L double{mustBeInteger};
    nn_level_L double{mustBeInteger} = -1;
    l_L double{mustBeFloat} = 0;
    m_L double{mustBeFloat} = 0;
    n_L double{mustBeFloat} = 0;
    options.overlap logical = false;
end
nhopping = length(L_1_L);
if nhopping == 1
    switch nargout
        case 2
            [varargout{1},varargout{2}] = ...
                HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
        case 1
            [varargout{1}] = ...
                HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
    end
    return;
end
exchange_list = L_1_L > L_2_L;
L_3_L = L_2_L(exchange_list);
L_2_L(exchange_list) = L_1_L(exchange_list);
L_1_L(exchange_list) = L_3_L;
m_3_L = m_1_L(exchange_list);
m_2_L(exchange_list) = m_1_L(exchange_list);
m_1_L(exchange_list) = m_3_L;
l_L(exchange_list) = - l_L(exchange_list);
m_L(exchange_list) = - m_L(exchange_list) ;
n_L(exchange_list) = - n_L(exchange_list);
Coeff = zeros(nhopping,3);
for i = 1:nhopping
    Coeff(i,:) = HR.TBSK_Coeff_gen(L_1_L(i),L_2_L(i),m_1_L(i),m_2_L(i),l_L(i),m_L(i),n_L(i));
end
nn_level_str_L= strcat(repmat('_',[nhopping,1]),string(nn_level_L));
if nargout == 1 && options.overlap
    varargout{1} = Coeff;
    return;
end
l2str_map=['s';'p';'d';'f'];
str_base_L = [l2str_map(L_1_L+1),l2str_map(L_2_L+1)];
if options.overlap
    Char0 = 'S';
else
    Char0 = 'V';
end
V_str_base_L = [repmat(Char0,[nhopping,1]),str_base_L];
V_str_1_L = strcat(V_str_base_L,repmat('S',[nhopping,1]),nn_level_str_L);
V_str_2_L = strcat(V_str_base_L,repmat('P',[nhopping,1]),nn_level_str_L);
V_str_3_L = strcat(V_str_base_L,repmat('D',[nhopping,1]),nn_level_str_L);

if nargout == 1
    E_hop_sym_out = sym(V_str_1_L,'real').*Coeff(:,1) +...
        sym(V_str_2_L,'real').*Coeff(:,2) +...
        sym(V_str_3_L,'real').*Coeff(:,3) ;
    varargout{1} = E_hop_sym_out;
    return;
else
    varargout{1} = Coeff;
    E_hop_sym_out = sym(V_str_1_L,'real').*Coeff(:,1) +...
        sym(V_str_2_L,'real').*Coeff(:,2) +...
        sym(V_str_3_L,'real').*Coeff(:,3) ;
    varargout{2} = E_hop_sym_out;
    return;
end


end
