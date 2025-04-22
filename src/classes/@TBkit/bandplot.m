function varargout = bandplot(TBkitobj,varargin)
% -------------- plot ------------------
if  TBkitobj.coe
    TBkitobj = TBkitobj.Subsall();
end
if  (nargin-1 ==1 && isvector(varargin{1}) ) 
    EIGENCAR = TBkitobj.EIGENCAR_gen();
    Ecut = varargin{1};
    varargin = varargin(2:end);
    bandplot(EIGENCAR,Ecut,TBkitobj.klist_l,TBkitobj.kpoints_l,TBkitobj.kpoints_name,varargin{:});
elseif nargin ==1 || ( (ischar(varargin{1}) || isstring(varargin{1}))  )
    EIGENCAR = TBkitobj.EIGENCAR_gen();
    Ecut=  [-3,3];
    bandplot(EIGENCAR,Ecut,TBkitobj.klist_l,TBkitobj.kpoints_l,TBkitobj.kpoints_name,varargin{:});
elseif ismatrix(varargin{1}) && ( (ischar(varargin{2}) || isstring(varargin{2}))  )
    EIGENCAR = varargin{1};
    varargin = varargin(2:end);
    Ecut=  [-3,3];
    varargout{:}  = bandplot(EIGENCAR,Ecut,TBkitobj.klist_l,TBkitobj.kpoints_l,TBkitobj.kpoints_name,varargin{:});
elseif  ismatrix(varargin{1}) && isnumeric(varargin{1}) && length(varargin{2}) ==2 
    EIGENCAR = varargin{1};Ecut = varargin{2};
    varargin = varargin(3:end);
    varargout{:} = bandplot(EIGENCAR,Ecut,TBkitobj.klist_l,TBkitobj.kpoints_l,TBkitobj.kpoints_name,varargin{:});
else
    varargout{:} = bandplot(varargin{:});
end
end