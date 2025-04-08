function H_hr = ForceToMat(H_hr)
switch H_hr.type
case 'sparse'
H_hr =H_hr.full();
case 'mat'
case 'list'
H_hr = H_hr.rewind;
otherwise
error('Not support yet.');
end
end
