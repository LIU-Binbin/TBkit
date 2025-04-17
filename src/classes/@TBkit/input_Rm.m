function TBkitobj = input_Rm(TBkitobj,Rm)
if nargin <2 && exist('POSCAR','file')
    [Rm,~,~,~,~] = POSCAR_read('POSCAR','vasp');
elseif nargin <2 && ~exist('POSCAR','file')
    Rm = eye(TBkitobj.Dim);
    %warning('POSCAR or Rm needed');
else
    if isa(Rm,'string') || isa(Rm,'char')
        [Rm,~,~,~,~] = POSCAR_read(Rm,'vasp');
    end
end
Rm = (Rm);
TBkitobj.Rm = Rm;
end