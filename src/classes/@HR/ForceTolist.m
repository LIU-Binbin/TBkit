function H_hr = ForceTolist(H_hr)
switch H_hr.type
case 'sparse'
H_hr =H_hr.full();
H_hr = H_hr.rewrite;
case 'mat'
H_hr = H_hr.rewrite;
case 'list'
otherwise
error('Not support yet.');
end
end
