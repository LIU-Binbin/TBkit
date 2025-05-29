function [H_soc_sym, lambda_syms] = soc_term_udud_add(elementL,qnumL,options)
arguments
    elementL = [];
    qnumL = [];
    options.mode {mustBeMember(options.mode, ["direct", "basis"])} = "direct";
    options.element_names = ["Mn", "Pt"];
    options.element_atom_nums = [2 2];
    options.element_projs = {2, [0,1,2]};
end

if strcmp(options.mode,'direct')
    element_names = options.element_names ;
    element_atom_nums = options.element_atom_nums ;
    element_projs = options.element_projs ;
    H_soc_p_uudd =  0.5 * ...
        [0 0 0 0 -1 1i;
        0 0 -1i 1 0 0;
        0 1i 0 -1i 0 0;
        0 1 1i 0 0 0;
        -1 0 0 0 0 1i;
        -1i 0 0 0 -1i 0];
    H_soc_d_uudd = 0.5 * ...
        [0 , 0 , 0 , 0 , 0 , 0 , -sqrt(3) , 1i*sqrt(3) , 0 , 0 ;
        0 , 0 , -1i , 0 , 0 , sqrt(3) , 0 , 0 , -1 , 1i ;
        0 , 1i , 0 , 0 , 0 , -1i*sqrt(3) , 0 , 0 , -1i , -1 ;
        0 , 0 , 0 , 0 , -2i , 0 , 1 , 1i , 0 , 0 ;
        0 , 0 , 0 , 2i , 0 , 0 , -1i , 1 , 0 , 0 ;
        0 , sqrt(3) , 1i*sqrt(3) , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
        -sqrt(3) , 0 , 0 , 1 , 1i , 0 , 0 , 1i , 0 , 0 ;
        -1i*sqrt(3) , 0 , 0 , -1i , 1 , 0 , -1i , 0 , 0 , 0 ;
        0 , -1 , 1i , 0 , 0 , 0 , 0 , 0 , 0 , 2i ;
        0 , -1i , -1 , 0 , 0 , 0 , 0 , 0 , -2i , 0];

    H_soc_s_udud = [0 0; 0 0]; % s-orbit have no SOC
    H_soc_p_udud = H_soc_p_uudd([1,4,2,5,3,6],[1,4,2,5,3,6]);
    H_soc_d_udud = H_soc_d_uudd([1,6,2,7,3,8,4,9,5,10],[1,6,2,7,3,8,4,9,5,10]);
    %%
    H_soc_sym = sym([]);
    for ie = 1:length(element_names)
        projs = element_projs{ie};
        lambda_ie = sym([]);

        for ia = 1:element_atom_nums(ie)
            for ip = 1:length(projs)
                switch projs(ip)
                    case 0
                        lambda = 0;
                        H_soc_sym = blkdiag(H_soc_sym, H_soc_s_udud * lambda);
                    case 1
                        lambda = sym("lambda__"+element_names(ie)+"_p");
                        % lambda = 0;
                        H_soc_sym = blkdiag(H_soc_sym, H_soc_p_udud * lambda);
                    case 2
                        lambda = sym("lambda__"+element_names(ie)+"_d");
                        H_soc_sym = blkdiag(H_soc_sym, H_soc_d_udud * lambda);
                    case 3
                        error("f orbit not supported yet")
                end
            end
        end
    end
elseif strcmp(options.mode,'basis')
    elements = readtable('elements.txt');
    EqnumL = [elementL,qnumL];
    [EqnumLsort,sort_label] = sortrows(EqnumL,[5],"descend");
    M  = TBkit.P2M(sort_label); 
    EQup = EqnumLsort(EqnumLsort(:,end) > 0,:);
    % L = [];
    for i = 1:size(EQup,1)
        ielement = EQup(i,1);
        sym_str = "lambda__" + elements(ielement,:).atom_symbol{1} + "_"+l2lname(EQup(i,3));
        sVar = str2sym(sym_str);
        %assume(symvar,'real');
        L(i,:) = Y_lm(EQup(i,3),EQup(i,4),sqrt(sVar),'An',ielement);
    end
    s = Spin(1/2);
    Hsoc = simplify( ...
        kron(s.Sz , L.Lz) + 1/2*( ...
        kron(s.Splus , L.Lminus) + ...
        kron(s.Sminus ,L.Lplus)));
    H_soc_sym = M'*Hsoc/M';
end

lambda_syms = symvar(H_soc_sym);
end

function name = l2lname(l)
switch l
    case 0
         name = "s";
    case 1
         name = "p";
    case 2   
        name = "d";
    case 3 
        name = "f";
    case -1   
        name = "sp";
    case -2
end
end



