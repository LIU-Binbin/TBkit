function OberserveValue = Observecar_gen(WAVECAR,Oper)
OberserveValue = zeros(size(WAVECAR,2),1);
for i = 1:size(WAVECAR,2)
    OberserveValue(i) = WAVECAR(:,i)'*Oper*WAVECAR(:,i);
end
end