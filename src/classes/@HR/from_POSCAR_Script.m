function H_hr = from_POSCAR_Script(varargin)
% FROM_POSCAR_SCRIPT Wrapper for POSCAR to HR conversion
%
%   H_hr = FROM_POSCAR_SCRIPT(varargin) acts as a wrapper for
%   HR.from_POSCAR_SE with variable input arguments.
%
%   INPUT ARGUMENTS:
%       varargin - Variable input arguments passed to from_POSCAR_SE
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object generated from POSCAR file
%
%   SEE ALSO:
%       HR.from_POSCAR_SE
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

H_hr = HR.from_POSCAR_SE(varargin{:});
end
