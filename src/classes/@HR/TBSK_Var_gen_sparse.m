function varargout = TBSK_Var_gen_sparse(L_1_L,L_2_L,m_1_L,m_2_L,Rlength_L,l_L,m_L,n_L,options)
arguments
L_1_L double{mustBeInteger};
L_2_L double{mustBeInteger};
m_1_L double{mustBeInteger};
m_2_L double{mustBeInteger};
Rlength_L double{mustBeNonnegative};
l_L  = 0;
m_L  = 0;
n_L  = 0;
options.overlap logical = false;
options.deltarule  double{mustBeMember(options.deltarule,[0,1,2])}= 0;
options.alpharule  double{mustBeMember(options.alpharule,[0,1,2])}= 0;
options.para = struct();
options.Rd double = -1;
end
nhopping = length(L_1_L);
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
l2str_map=['s';'p';'d';'f'];
str_base_L = [l2str_map(L_1_L+1),l2str_map(L_2_L+1)];
V_str_base_L = [repmat('V',[nhopping,1]),str_base_L];
V_str_1_L = string(strcat(V_str_base_L,repmat('S',[nhopping,1])));
V_str_2_L = string(strcat(V_str_base_L,repmat('P',[nhopping,1])));
V_str_3_L = string(strcat(V_str_base_L,repmat('D',[nhopping,1])));
strvar = options.para.strvar;
numvar = options.para.numvar;
[V_seq_1_L_contain,V_seq_1_L] = ismember(V_str_1_L,strvar);
[V_seq_2_L_contain,V_seq_2_L] = ismember(V_str_2_L,strvar);
[V_seq_3_L_contain,V_seq_3_L] = ismember(V_str_3_L,strvar);
V_seq_1_L = V_seq_1_L(V_seq_1_L_contain);
V_seq_2_L = V_seq_2_L(V_seq_2_L_contain);
V_seq_3_L = V_seq_3_L(V_seq_3_L_contain);
V_num_1_L = zeros(nhopping,1);V_num_2_L = V_num_1_L;V_num_3_L = V_num_1_L;
V_num_1_L(V_seq_1_L_contain) = numvar(V_seq_1_L);
V_num_2_L(V_seq_2_L_contain) = numvar(V_seq_2_L);
V_num_3_L(V_seq_3_L_contain) = numvar(V_seq_3_L);
V_scale1_L = ones(nhopping,1);V_scale2_L = V_scale1_L;V_scale3_L = V_scale1_L;
if options.alpharule > 0
if options.alpharule == 1 && length(options.para.delta) == 1
V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.Rd(V_seq_1_L)*options.para.delta);
V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.Rd(V_seq_2_L)*options.para.delta);
V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.Rd(V_seq_3_L)*options.para.delta);
elseif options.alpharule == 2 && length(options.para.delta) > 1
V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.Rd(V_seq_1_L).*options.para.delta(V_seq_1_L).');
V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.Rd(V_seq_2_L).*options.para.delta(V_seq_2_L).');
V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.Rd(V_seq_3_L).*options.para.delta(V_seq_3_L).');
end
end
if options.deltarule > 0
if options.deltarule == 1 && length(options.para.delta) == 1
V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')/options.para.delta);
V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')/options.para.delta);
V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')/options.para.delta);
elseif options.deltarule == 2 && length(options.para.delta) > 1
V_scale1_L(V_seq_1_L_contain) = exp(-(Rlength_L(V_seq_1_L_contain) - options.Rd(V_seq_1_L).')./options.para.delta(V_seq_1_L).');
V_scale2_L(V_seq_2_L_contain) = exp(-(Rlength_L(V_seq_2_L_contain) - options.Rd(V_seq_2_L).')./options.para.delta(V_seq_2_L).');
V_scale3_L(V_seq_3_L_contain) = exp(-(Rlength_L(V_seq_3_L_contain) - options.Rd(V_seq_3_L).')./options.para.delta(V_seq_3_L).');
end
end
HoppingL = Coeff(:,1).*V_num_1_L.*V_scale1_L + ...
Coeff(:,2).*V_num_2_L.*V_scale2_L + ...
Coeff(:,3).*V_num_3_L.*V_scale3_L;
if nargout ==1
varargout{1} = HoppingL;
end
end
