function [H_soc_sym, lambda_syms] = soc_term_udud_add(element_names, element_atom_nums, element_projs)
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
H_soc_sym = [];
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
                    % lambda = sym("lambda_"+element_names(ie)+"_p");
                    lambda = 0;
                    H_soc_sym = blkdiag(H_soc_sym, H_soc_p_udud * lambda);
                case 2
                    lambda = sym("lambda_"+element_names(ie)+"_d");
                    H_soc_sym = blkdiag(H_soc_sym, H_soc_d_udud * lambda);
                case 3
                    error("f orbit not supported yet")
            end
        end
    end
end
lambda_syms = symvar(H_soc_sym);
end



