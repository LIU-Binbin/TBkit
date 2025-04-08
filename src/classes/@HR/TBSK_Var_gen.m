function varargout = TBSK_Var_gen(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,options)
arguments
L_1_L double{mustBeInteger};
L_2_L double{mustBeInteger};
m_1_L double{mustBeInteger};
m_2_L double{mustBeInteger};
nn_level_L double{mustBeInteger} = -1;
l_L  = 0;
m_L  = 0;
n_L  = 0;
options.overlap logical = false;
end
nhopping = length(L_1_L);
if nhopping == 1
switch nargout
case 3
[varargout{1},varargout{2},varargout{3}] = ...
HR.TBSK_Var_gen_single(L_1_L,L_2_L,m_1_L,m_2_L,nn_level_L,l_L,m_L,n_L,'overlap', options.overlap );
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
V_str_base_L = [repmat('V',[nhopping,1]),str_base_L];
V_str_1_L = strcat(V_str_base_L,repmat('S',[nhopping,1]),nn_level_str_L);
V_str_2_L = strcat(V_str_base_L,repmat('P',[nhopping,1]),nn_level_str_L);
V_str_3_L = strcat(V_str_base_L,repmat('D',[nhopping,1]),nn_level_str_L);
if options.overlap
S_str_base_L = [repmat('S',[nhopping,1]),str_base_L];
S_str_1_L = strcat(S_str_base_L,repmat('S',[nhopping,1]),nn_level_str_L);
S_str_2_L = strcat(S_str_base_L,repmat('P',[nhopping,1]),nn_level_str_L);
S_str_3_L = strcat(S_str_base_L,repmat('D',[nhopping,1]),nn_level_str_L);
end
if nargout == 3
varargout{1} = Coeff;
varargout{2} = [V_str_1_L,V_str_2_L,V_str_3_L];
varargout{3} = [S_str_1_L,S_str_2_L,S_str_3_L];
return;
end
E_hop_sym_out = sym(V_str_1_L,'real').*Coeff(:,1) +...
sym(V_str_2_L,'real').*Coeff(:,2) +...
sym(V_str_3_L,'real').*Coeff(:,3) ;
if options.overlap
E_integral_sym_out = sym(S_str_1_L,'real').*Coeff(:,1) +...
sym(S_str_2_L,'real').*Coeff(:,2) +...
sym(S_str_3_L,'real').*Coeff(:,3) ;
end
switch nargout
case 1
varargout{1} = E_hop_sym_out;
case 2
varargout{1} = E_hop_sym_out;
varargout{2} = E_integral_sym_out;
otherwise
error('!!');
end
end
