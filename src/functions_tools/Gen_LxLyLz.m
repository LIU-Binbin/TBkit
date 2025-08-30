function Gen_LxLyLz()
[orbL, elementL, quantumL] = wout_read();
%%
ZeroMat = zeros(size(orbL,1));
for i = 1:size(orbL,1)
    for j = 1:size(orbL,1)
        ZeroMat(i,j) = ~all(orbL(i,:)==orbL(j,:));
    end
end

EqnumL = [elementL,quantumL];
[EqnumLsort,sort_label] = sortrows(EqnumL,5,"descend");
M  = TBkit.P2M(sort_label); 
EQup = EqnumLsort(EqnumLsort(:,end) > 0,:);

clear L
for i = 1:size(EQup,1)
    ielement = EQup(i,1);
    %assume(symvar,'real');
    L(i,:) = Y_lm(EQup(i,3), EQup(i,4) ,'An',ielement);
end

I2 = sym(eye(2));

Lx = kron(I2, (L.Lplus+L.Lminus)/sym(2));
Ly = kron(I2, (L.Lplus-L.Lminus)/sym(2i));
Lz = kron(I2, (L.Lz));

Lx = M'*Lx/M';
Ly = M'*Ly/M';
Lz = M'*Lz/M';
% 
Lx(logical(ZeroMat)) = 0;
Ly(logical(ZeroMat)) = 0;
Lz(logical(ZeroMat)) = 0;
%% wannier tools use inner uudd basis
WAN_NUM = length(Lx);

udud2uudd = [1:2:(WAN_NUM-1),2:2:(WAN_NUM)];
Lx = Lx(udud2uudd,udud2uudd);
Ly = Ly(udud2uudd,udud2uudd);
Lz = Lz(udud2uudd,udud2uudd);

%% write to file
fid1 = fopen('Lx.dat', 'w');
fid2 = fopen('Ly.dat', 'w');
fid3 = fopen('Lz.dat', 'w');
for k = 1:numel(Lx)
    fprintf(fid1, '%.16e %.16e\n', real(Lx(k)), imag(Lx(k)));
    fprintf(fid2, '%.16e %.16e\n', real(Ly(k)), imag(Ly(k)));
    fprintf(fid3, '%.16e %.16e\n', real(Lz(k)), imag(Lz(k)));
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
end
