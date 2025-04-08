function H_hr = cut_orb(H_hr,rm_list,options)
arguments
H_hr HR;
rm_list = [];
options.rmfunc  function_handle=@()(1);
end
orb_tmp = H_hr.orbL;
if isempty(rm_list) && ~strcmp(functions(options.rmfunc).function , '@()(1)')
rm_list = options.rmfunc(orb_tmp(:,1),orb_tmp(:,2),orb_tmp(:,3));
elseif isempty(rm_list)
rm_list = zeros(size(H_hr.orbL,1),1);
else
rmlist = zeros(size(H_hr.orbL,1),1);
rmlist(rm_list) = 1;
rm_list = rmlist;
end
wan_list = ~rm_list;
H_hr = H_hr.reseq(wan_list);
end
