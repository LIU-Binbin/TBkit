function Gen_LS(hr)
arguments
    hr = []
end
if isempty(hr)
    [orbL, elementL, quantumL] = wout_read();
else
    orbL = hr.orbL;
    elementL = hr.elementL;
    quantumL = hr.quantumL;
end
%%
[H_soc_sym, lambda_syms] = SOC_on_site(orbL, elementL, quantumL);
LS_mat_num = double(subs(H_soc_sym, lambda_syms, ones(1,length(lambda_syms))));
%% wannier tools use inner uudd basis
WAN_NUM = length(LS_mat_num);

udud2uudd = [1:2:(WAN_NUM-1),2:2:(WAN_NUM)];
LS_mat_num = LS_mat_num(udud2uudd,udud2uudd);
%%
fid1 = fopen('LS.dat', 'w');
for k = 1:numel(LS_mat_num)
    fprintf(fid1, '%.16e %.16e\n', real(LS_mat_num(k)), imag(LS_mat_num(k)));
end
fclose(fid1);
end