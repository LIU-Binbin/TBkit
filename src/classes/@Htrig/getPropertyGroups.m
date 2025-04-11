function propgrp = getPropertyGroups(~)
% getPropertyGroups Creates a PropertyGroup object containing core properties of the Hamiltonian system.
% This function defines a set of grouped properties for the Hamiltonian structure,
% typically used within a class to organize related attributes for tooling purposes.
%
% propgrp = getPropertyGroups(~) returns a PropertyGroup object that groups the 
% following properties: Basis_num, Kinds, Type, HsymL, symvar_list, and Dim. These 
% properties are essential for defining and manipulating the Hamiltonian's 
% configuration within MATLAB's class framework.
%
% Inputs:
%   ~: Unused input parameter (for compatibility with parent class requirements).
%
% Outputs:
%   propgrp: PropertyGroup object containing the specified properties. This object
%            is used by MATLAB's class property tools to group related attributes.
%
% Example:
%   % Access grouped properties in a class display
%   propgrp = getPropertyGroups(~);
%   showProperties(propgrp);

    proplist = {'Basis_num','Kinds','Type','HsymL','symvar_list','Dim'};
    propgrp = matlab.mixin.util.PropertyGroup(proplist);
end
