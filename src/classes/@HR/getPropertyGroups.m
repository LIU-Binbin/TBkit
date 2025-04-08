function propgrp = getPropertyGroups(~)
proplist = {'WAN_NUM','NRPTS','Type','HcoeL','HnumL','vectorL'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
