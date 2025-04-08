function H_hr = set_hop_single(H_hr,amp,hi,hj,vector,mode)
V = H_hr.vectorL;
if strcmp(H_hr.Type,'list')
vector = [vector,hi,hj];
else
end
if isempty(V)
seq = H_hr.NRPTS +1;
H_hr = H_hr.add_empty_one(vector);
else
[~,seq]=ismember((vector),V,'rows');
if seq == 0
seq = H_hr.NRPTS +1;
H_hr = H_hr.add_empty_one(vector);
end
end
switch H_hr.Type
case 'list'
switch mode
case 'set'
if H_hr.HnumL(seq) ~= 0
warning('May be you should use add mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell(vector(1,:));
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HnumL(seq) = amp ;
case 'add'
H_hr.HnumL(seq) = H_hr.HnumL(seq) + amp  ;
case 'sym'
if H_hr.HcoeL(seq) ~= sym(0)
warning('May be you should use symadd mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell(vector(1,:));
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HcoeL(seq) =  amp  ;
case 'symadd'
H_hr.HcoeL(seq) = H_hr.HcoeL(seq) + amp  ;
end
case 'sparse'
switch mode
case 'set'
if H_hr.HnumL{seq}(hi,hj) ~= 0
warning('May be you should use add mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell([vector(1,:),hi,hj]);
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HnumL{seq}(hi,hj) = amp ;
case 'add'
H_hr.HnumL{seq}(hi,hj) = H_hr.HnumL(hi,hj,seq) + amp  ;
case 'sym'
if H_hr.HcoeL{seq}(hi,hj) ~= sym(0)
warning('May be you should use symadd mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell([vector(1,:),hi,hj]);
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HcoeL{seq}(hi,hj) =  amp  ;
case 'symadd'
H_hr.HcoeL{seq}(hi,hj) = H_hr.HcoeL{seq}(hi,hj) + amp  ;
end
otherwise
switch mode
case 'set'
if H_hr.HnumL(hi,hj,seq) ~= 0
warning('May be you should use add mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell([vector(1,:),hi,hj]);
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HnumL(hi,hj,seq) = amp ;
case 'add'
H_hr.HnumL(hi,hj,seq) = H_hr.HnumL(hi,hj,seq) + amp  ;
case 'sym'
if H_hr.HcoeL(hi,hj,seq) ~= sym(0)
warning('May be you should use symadd mode on this NRPT and hi hj');
Format = fold(@strcat,repmat("
Input = num2cell([vector(1,:),hi,hj]);
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HcoeL(hi,hj,seq) =  amp  ;
case 'symadd'
H_hr.HcoeL(hi,hj,seq) = H_hr.HcoeL(hi,hj,seq) + amp  ;
end
end
end
