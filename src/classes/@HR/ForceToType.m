function H_hr = ForceToType(H_hr,Type)
switch Type
case 'sparse'
H_hr = H_hr.ForceTosparse();
case 'mat'
H_hr = H_hr.ForceToMat();
case 'list'
H_hr = H_hr.ForceTolist();
otherwise
error('Not support yet.');
end
end
