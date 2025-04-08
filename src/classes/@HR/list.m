function A = list(H_hr,vectorL,options)
arguments
H_hr;
vectorL = [0,0,0];
options.vpa =true;
options.digits = 5;
options.numeric = false;
options.disp = true;
end
if H_hr.vectorhopping
H_hr = H_hr.GenfromOrth();
end
H_hr = H_hr.rewrite();
RFORMAT =  fold(@strcat,"R"+string([1:H_hr.Dim])+" ");
FORMAT = strcat('# ',RFORMAT,'i j real imag\n');
fprintf(FORMAT);
if nargin <2
if H_hr.coe && ~options.numeric
if options.vpa
A = [sym(H_hr.vectorL),...
vpa(real(H_hr.HcoeL), options.digits) ,...
vpa(imag(H_hr.HcoeL),options.digits)];
else
A = [sym(H_hr.vectorL),...
(real(H_hr.HcoeL)) ,...
imag(H_hr.HcoeL)];
end
disp(A);
elseif H_hr.num
disp([double(H_hr.vectorL),real(H_hr.HnumL),imag(H_hr.HnumL)]);
end
else
vectorList = double(H_hr.vectorL);
if isa(vectorL,'double')
switch size(vectorL,2)
case 2
[seq] = find(all(vectorL == vectorList(:,H_hr.Dim+1:H_hr.Dim+2),2));
case 3
[seq] = find(all(vectorL == vectorList(:,1:H_hr.Dim),2));
case 5
[seq] = find(all(vectorL == vectorList,2));
end
elseif isa(vectorL,'sym')
[seq] = find(park.strcontain(string(H_hr.HcoeL),string(vectorL)));
end
if H_hr.coe  && ~options.numeric
if options.vpa
A = [sym(H_hr.vectorL(seq,:)),...
vpa(real(H_hr.HcoeL(seq,:)), options.digits) ,...
vpa(imag(H_hr.HcoeL(seq,:)),options.digits)];
else
A = [sym(H_hr.vectorL(seq,:)),real(H_hr.HcoeL(seq,:)) ,imag(H_hr.HcoeL(seq,:))];
end
if options.disp
disp(vpa(A));
end
elseif H_hr.num
disp([double(H_hr.vectorL(seq,:)),real(H_hr.HnumL(seq,:)),imag(H_hr.HnumL(seq,:))]);
end
end
end
