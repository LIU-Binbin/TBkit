function propgrp = getPropertyGroups(~)
proplist = {'Basis_num','Kinds','Type','HsymL','symvar_list','Dim'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
