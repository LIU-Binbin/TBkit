function propgrp = getPropertyGroups(~)
%GETPROPERTYGROUPS Custom property display for HK class
%
% Syntax:
%   propgrp = getPropertyGroups(~)
%
% Output:
%   propgrp - PropertyGroup object controlling display
%
% Description:
%   Defines which properties are displayed when viewing an HK object.
%   Controls MATLAB's default property grouping behavior by specifying
%   a subset of key properties to show (Degree, Kinds, HsymL, HcoeL, HnumL).
%
% Note:
%   Part of MATLAB's CustomDisplay interface - not meant for direct calls
%
% Implementation:
%   Uses matlab.mixin.util.PropertyGroup to create display group
proplist = {'Degree','Kinds','HsymL','HcoeL','HnumL'};
propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
