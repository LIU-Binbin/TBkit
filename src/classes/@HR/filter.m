function H_hr = filter(H_hr,Accuracy)
if nargin < 2
Accuracy = 1e-6;
end
switch H_hr.Type
case {'list','mat'}
HnumList = H_hr.HnumL;
HnumList_real = real(HnumList);
HnumList_imag = imag(HnumList);
HnumList_real(abs(HnumList_real) < Accuracy) = 0;
HnumList_imag(abs(HnumList_imag) < Accuracy) = 0;
H_hr.HnumL = HnumList_real + 1i*HnumList_imag;
otherwise
end
end
