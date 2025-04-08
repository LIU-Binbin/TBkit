function H_hr = input_orb_struct(H_hr,filename,mode,options)
arguments
H_hr HR;
filename string ='POSCAR';
mode char {mustBeMember(mode,{'vasp','tbsk','sym'})} = 'vasp';
options.symbolic logical = false;
options.Operation logical = false;
options.warning logical = true
options.spin  char {mustBeMember(options.spin,{'spinless','wannier','block'})}= 'spinless';
end
optionsCell = namedargs2cell(options);
H_hr.Basis_num = H_hr.WAN_NUM;
H_hr = input_orb_struct@TBkit(H_hr,filename,mode,...
optionsCell{:});
end
