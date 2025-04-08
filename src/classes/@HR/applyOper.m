function H_hr = applyOper(H_hr,SymOper,options)
arguments
H_hr HR;
SymOper Oper = Oper();
options.generator = false;
options.HcoeL = false;
options.Accuracy = 1e-6;
options.fast = false;
options.center = [0,0,0];
options.Ugen = false;
end
options2 =options;
nSymOper = length(SymOper);
if options.Ugen
try
BasisFunction = BasisFunc(H_hr);
SymOper = SymOper.Ugen(BasisFunction,'Rm',H_hr.Rm,'center',options.center);
catch
end
end
if options.fast
if isempty(H_hr.AvectorL) && isempty(H_hr.BvectorL)
H_hr = H_hr.init('fast',true);
H_hr = H_hr.hermitize();
end
if options.generator
options2.generator = false;
optionsCell = namedargs2cell(options2);
for i = 1:nSymOper
fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
disp(SymOper(i));
SymOper_tmp = SymOper(i).generate_group();
H_hr =  applyOper(H_hr,SymOper_tmp,optionsCell{:});
fprintf('----------   SymVarNum: %d   ----------\n',...
rank(H_hr.CvectorL));
end
H_hr = H_hr.hermitize;
H_hr = H_hr.simplify(options.Accuracy);
else
H_hr_R = H_hr;
nSymOper = length(SymOper);
pb = TBkit_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
for j = 1:nSymOper
[H_hr_R(j),H_hr] = applyRU(H_hr,SymOper(j));
pb.print(j,nSymOper);
end
pb.delete();
H_hr = sum(H_hr_R);
H_hr = H_hr.hermitize;
H_hr = H_hr.simplify(options.Accuracy);
end
return;
end
if ~H_hr.coe && ~H_hr.num
H_hr = H_hr.init();
H_hr = H_hr.hermitize();
end
if ~strcmp(H_hr.Type , 'list')
H_hr = H_hr.rewrite();
end
if length(SymOper) == 1
if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
return;
end
nSymOper = length(SymOper);
fprintf('******** apply (%d/%d)symmetry ********\n',1,nSymOper);
disp(SymOper);
if options.generator
SymOper_tmp = SymOper.generate_group();
nSymOper_tmp = length(SymOper_tmp);
pb = TBkit_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
H_hr_R = H_hr;
for j = 1:nSymOper_tmp
pb.print(j,nSymOper_tmp);
[H_hr_R(j),H_hr] = applyRU(H_hr,SymOper_tmp(j));
end
pb.delete();
H_hr = sum(H_hr_R)/nSymOper_tmp;
H_hr = H_hr.simplify();
elseif options.HcoeL
[H_hr_R,H_hr] = applyRU(H_hr,SymOper);
if isequal(H_hr_R.HcoeL,H_hr.HcoeL)
else
Equationlist_r = (real(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
Equationlist_i = (imag(H_hr.HcoeL - H_hr_R.HcoeL) == 0);
Equationlist_r = HR.isolateAll(Equationlist_r);
Equationlist_i = HR.isolateAll(Equationlist_i);
HcoeLtmp = H_hr.HcoeL ;
HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
H_hr.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
end
H_hr = H_hr.simplify();
end
else
nSymOper = length(SymOper);
for i = 1:nSymOper
fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
disp(SymOper(i));
if options.generator
SymOper_tmp = SymOper(i).generate_group();
nSymOper_tmp = length(SymOper_tmp);
pb = TBkit_tool_outer.CmdLineProgressBar('Applying Symmetry ...');
H_hr_R = H_hr;
for j = 1:nSymOper_tmp
pb.print(j,nSymOper_tmp);
[H_hr_R(j),H_hr] = applyRU(H_hr,SymOper_tmp(j));
end
pb.delete();
fprintf('Debug');
H_hr = sum(H_hr_R)/nSymOper_tmp;
H_hr = H_hr.simplify(options.Accuracy);
else
%fprintf('    ');
H_hr = H_hr.applyOper(SymOper(i),'generator','false');
end
fprintf('----------   SymVarNum: %d   ----------\n',length(H_hr.symvar_list));
end
end
end
