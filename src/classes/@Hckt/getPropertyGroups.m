function propgrp = getPropertyGroups(~)
proplist = {'title','dim','Nports','Components','vectorAll','Netlist'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
