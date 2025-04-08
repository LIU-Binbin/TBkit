function H_hr = descritize(H_hr,Nslab,options)
arguments
H_hr HR;
Nslab double = [0,10,0];
options.rmfunc function_handle=@()(1);
options.Rotation  = sym(eye(3));
options.np = 1;
options.vacuum_mode = 0;
end
if strcmp(functions(options.rmfunc).function , '@()(1)')
rm_mode =false;
else
rm_mode =true;
end
if isequal(options.Rotation ,sym(eye(3)))
rotate_mode = false;
else
rotate_mode = true;
end
end
