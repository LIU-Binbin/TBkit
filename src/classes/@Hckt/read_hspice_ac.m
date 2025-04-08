function [simulation_result]=read_hspice_ac(filename,options)
arguments
filename = "Hckt.ac0";
options.fast = true;
options.filenameInformation = '';
options.ErrorMode = 0;
end
optionscell = namedargs2cell(options);
[simulation_result]=Hckt.read_hspice_tr(filename,optionscell{:});
end
