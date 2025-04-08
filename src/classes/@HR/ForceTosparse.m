function H_hr = ForceTosparse(H_hr)
switch H_hr.type
case 'sparse'
case 'mat'
H_hr = H_hr.sparse;
case 'list'
H_hr = H_hr.rewind;
H_hr = H_hr.rewind;
otherwise
error('Not support yet.');
end
end
