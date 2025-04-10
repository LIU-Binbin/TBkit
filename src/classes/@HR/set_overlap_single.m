function H_hr = set_overlap_single(H_hr,amp,si,sj,vector,mode)
% SET_OVERLAP_SINGLE Set single overlap element in HR object
%
%   H_HR = SET_OVERLAP_SINGLE(H_HR,AMP,SI,SJ,VECTOR,MODE) sets single overlap
%
%   Inputs:
%       H_hr - HR object to modify
%       amp - Overlap amplitude
%       si - Initial orbital index
%       sj - Final orbital index
%       vector - Lattice vector for overlap
%       mode - Operation mode ('set', 'add', 'sym', 'symadd')
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - Low-level function for overlap terms
%       - Similar to set_hop_single but for overlaps
%       - Provides warnings for potential overwrites
V = H_hr.vectorL_overlap;
if strcmp(H_hr.Type,'list')
vector = [vector,si,sj];
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
if H_hr.SnumL(seq) ~= 0
warning('May be you should use add mode on this NRPT and hi hj');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(H_hr.Dim+2));
end
H_hr.SnumL(seq) = amp ;
case 'add'
H_hr.SnumL(seq) = H_hr.SnumL(seq) + amp  ;
case 'sym'
if H_hr.ScoeL(seq) ~= sym(0)
warning('May be you should use symadd mode on this NRPT and hi hj');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),vector(4),vector(H_hr.Dim+2));
end
H_hr.ScoeL(seq) =  amp  ;
case 'symadd'
H_hr.ScoeL(seq) = H_hr.ScoeL(seq) + amp  ;
end
otherwise
switch mode
case 'set'
if H_hr.SnumL(si,sj,seq) ~= 0
warning('May be you should use add mode on this NRPT and si sj');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),si,sj);
end
H_hr.SnumL(si,sj,seq) = amp ;
case 'add'
H_hr.SnumL(si,sj,seq) = H_hr.SnumL(si,sj,seq) + amp  ;
case 'sym'
if H_hr.ScoeL(si,sj,seq) ~= sym(0)
warning('May be you should use symadd mode on this NRPT and si sj');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),si,sj);
end
H_hr.ScoeL(si,sj,seq) = amp  ;
case 'symadd'
H_hr.ScoeL(si,sj,seq) = H_hr.ScoeL(si,sj,seq) + amp  ;
end
end
end
