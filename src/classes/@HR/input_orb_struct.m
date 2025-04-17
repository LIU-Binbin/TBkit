function H_hr = input_orb_struct(H_hr,filename,mode,options)
%INPUT_ORB_STRUCT Import orbital structure into HR object
%   Reads orbital structure from file and imports it into an HR object.
%
%   Syntax:
%       H_hr = input_orb_struct(H_hr)
%       H_hr = input_orb_struct(H_hr,filename)
%       H_hr = input_orb_struct(H_hr,filename,mode)
%       H_hr = input_orb_struct(H_hr,filename,mode,options)
%
%   Inputs:
%       H_hr - Target HR object
%       filename - Input filename (default: 'POSCAR')
%       mode - Input mode ('vasp', 'tbsk', or 'sym') (default: 'vasp')
%       options - Optional name-value pairs:
%           symbolic - Logical flag for symbolic mode (default: false)
%           Operation - Logical flag for operations (default: false)
%           warning - Logical flag for warnings (default: true)
%           spin - Spin configuration ('spinless','wannier','block')
%
%   Outputs:
%       H_hr - HR object with imported orbital structure
%
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
H_hr = input_orb_struct@TBkit(H_hr,filename,mode,optionsCell{:});
end
