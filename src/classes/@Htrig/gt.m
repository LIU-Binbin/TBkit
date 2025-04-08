function C = gt(B,A)
if isa(A,'Htrig') && isa(B,'Htrig')
H_htrig1 = A;
H_htrig2 = B;
if H_htrig1.Bassis_num ~= H_htrig2.Bassis_num
error('Bassis_num different');
end
error('not support at present.')
elseif isa(A,'Htrig') && ~isa(B,'Htrig')
switch class(B)
case 'char'
switch B(1)
case {'P','p'}
C = A.input_orb_struct(B,'sym');
C.Rm = sym(C.Rm );
C.orbL = sym(C.orbL );
case {'w','W'}
case {'k','K'}
C = A.kpathgen3D(B);
otherwise
end
case 'double'
switch size(B,1)
case A.WAN_NUN
C = A.input_orb_init(B);
otherwise
end
otherwise
end
elseif ~isa(A,'Htrig') && ~isa(B,'Htrig')
error('not support at present.');
end
end
