% bandplot.m
function varargout = bandplot(TBkitobj, varargin)
%BANDPLOT Plot band structure from TBkit object
%   BANDPLOT(TBkitobj, ...) generates band structure plots from a TBkit
%   object. Supports multiple input formats for energy cutoff and plotting
%   parameters.
%
% Input Arguments:
%   TBkitobj : TBkit object
%       Tight-binding calculation object containing k-point data
%   varargin : Flexible arguments
%       - Ecut: Energy range [min, max] (default: [-3, 3])
%       - EIGENCAR: Precomputed eigenvalue matrix
%       - Plotting parameters (line style, color, etc.)
%
% Output Arguments:
%   varargout : Handle array
%       Graphics handles when using EIGENCAR input
%
% Overloaded Syntax Support:
%   1. bandplot(TBkitobj)
%   2. bandplot(TBkitobj, Ecut)
%   3. bandplot(TBkitobj, EIGENCAR)
%   4. bandplot(TBkitobj, EIGENCAR, Ecut)
%
% Example:
%   TBobj = TBkit(...);
%   bandplot(TBobj, [-2, 2], 'Color','red');

% Process symbolic coefficients if needed
if TBkitobj.coe
    TBkitobj = TBkitobj.Subsall();
end

% Handle different input configurations
if (nargin-1 == 1 && isvector(varargin{1}))
    % Case: bandplot(TBkitobj, Ecut)
    EIGENCAR = TBkitobj.EIGENCAR_gen();
    Ecut = varargin{1};
    varargin = varargin(2:end);
    bandplot(EIGENCAR, Ecut, TBkitobj.klist_l, TBkitobj.kpoints_l, TBkitobj.kpoints_name, varargin{:});
elseif nargin == 1 || ((ischar(varargin{1}) || isstring(varargin{1}))
    % Case: bandplot(TBkitobj) or bandplot(TBkitobj, style...)
    EIGENCAR = TBkitobj.EIGENCAR_gen();
    Ecut = [-3, 3];
    bandplot(EIGENCAR, Ecut, TBkitobj.klist_l, TBkitobj.kpoints_l, TBkitobj.kpoints_name, varargin{:});
elseif ismatrix(varargin{1}) && ((ischar(varargin{2}) || isstring(varargin{2}))
    % Case: bandplot(TBkitobj, EIGENCAR, style...)
    EIGENCAR = varargin{1};
    varargin = varargin(2:end);
    Ecut = [-3, 3];
    varargout{:} = bandplot(EIGENCAR, Ecut, TBkitobj.klist_l, TBkitobj.kpoints_l, TBkitobj.kpoints_name, varargin{:});
elseif ismatrix(varargin{1}) && isnumeric(varargin{1}) && length(varargin{2}) == 2
    % Case: bandplot(TBkitobj, EIGENCAR, Ecut)
    EIGENCAR = varargin{1};
    Ecut = varargin{2};
    varargin = varargin(3:end);
    varargout{:} = bandplot(EIGENCAR, Ecut, TBkitobj.klist_l, TBkitobj.kpoints_l, TBkitobj.kpoints_name, varargin{:});
else
    % Generic case forwarding
    varargout{:} = bandplot(varargin{:});
end
end