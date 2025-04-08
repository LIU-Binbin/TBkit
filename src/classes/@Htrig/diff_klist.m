function varargout = diff_klist(H_htrig,dir,klist,options)
arguments
H_htrig
dir = 1;
klist = [];
options.Accuracy = 1e-6;
end
if ~strcmp(H_htrig.Type,'list')
end
optionscell = namedargs2cell(options);
H_htrig_tmp = H_htrig;
varargout{nargout} = H_htrig_tmp.diff(dir(numel(dir)),optionscell{:}).HCAR_gen(klist);
for i = 1:numel(dir)-1
H_htrig_tmp = H_htrig;
varargout{i} = H_htrig_tmp.diff(dir(i),optionscell{:}).HCAR_gen(klist);
end
end
