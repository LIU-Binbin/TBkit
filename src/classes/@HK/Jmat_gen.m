function J_mat = Jmat_gen(H_hk,QuantumL,options)
arguments
H_hk  HK;
QuantumL double = H_hk.quantumL;
options.include_0 logical = false;
options.strict logical = false;
end
J_mat = Oper.spin_matrices_from_orb(QuantumL,options.include_0,options.strict);
end
