function propgrp = getPropertyGroups(~)
proplist = {'Degree','Kinds','HsymL','HcoeL','HnumL'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
