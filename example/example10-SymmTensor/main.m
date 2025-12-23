%[text] # Symmetry Adapted Tensor
%[text] ## Examples for specifying the operators manually
% P = Oper.inversion();
% Mx = Oper.mirror([1,0,0]);
% My = Oper.mirror([0,1,0]);
% Mz = Oper.mirror([0,0,1]);
% C6x = Oper.rotation(1/6,[1,0,0]);
% C6z = Oper.rotation(1/6,[0,0,1]);
% C3x = Oper.rotation(1/3,[1,0,0]);
% C3z = Oper.rotation(1/3,[0,0,1]);
% S6x = Mx*C6x;
% S6z = Mz*C6z;
%%
Msg = read_Magnetic_Sym_El([166,98]);
Gen_list = Msg;

% AHC_1st     = 'a{V2}';
% AHC_2nd_QMD = 'a{V2}V';
% AHC_2nd_BCD = '{V2}V';
% SHC         = 'eV3';

jahn_symbol_Str = 'eV3';
Tensor = jahn_symbol(jahn_symbol_Str);

for i = 1:length(Gen_list)
    Tensor = group_transformation(Tensor, Gen_list(i));
end

[~, ~, SymMatDisplay] = pretty(Tensor,"Table"); %[output:7ea1fcb1]
SymMatDisplay %[output:16c11ae5]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:7ea1fcb1]
%   data: {"dataType":"text","outputData":{"text":"Independent Elements: 4\n","truncated":false}}
%---
%[output:16c11ae5]
%   data: {"dataType":"symbolic","outputData":{"name":"SymMatDisplay","value":"\\left(\\begin{array}{cccccccccc}\n\\mathrm{Tensor} & 11 & 21 & 31 & 12 & 22 & 32 & 13 & 23 & 33\\\\\n1 & \\chi_{1,1,1}  & 0 & 0 & 0 & -\\chi_{1,1,1}  & -\\chi_{2,3,1}  & 0 & -\\chi_{2,1,3}  & 0\\\\\n2 & 0 & -\\chi_{1,1,1}  & \\chi_{2,3,1}  & -\\chi_{1,1,1}  & 0 & 0 & \\chi_{2,1,3}  & 0 & 0\\\\\n3 & 0 & \\chi_{3,2,1}  & 0 & -\\chi_{3,2,1}  & 0 & 0 & 0 & 0 & 0\n\\end{array}\\right)"}}
%---
