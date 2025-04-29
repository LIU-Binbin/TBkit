function TBkitobj = TBkitCopy(TBkitobj,TBkitobj_in)
 %TBKITCOPY Copy properties between TBkit objects
    %
    %   Syntax:
    %       TBkitobj = TBkitCopy(TBkitobj,TBkitobj_in)
    %
    %   Description:
    %       Copies specified properties from one TBkit object to another.
    %
    %   Inputs:
    %       TBkitobj    - Destination TBkit object
    %       TBkitobj_in - Source TBkit object
    %
    %   Output:
    %       TBkitobj - Modified destination TBkit object
    CopyItem = [...
    "Rm","orbL","Dim","elementL","quantumL","orb_symL","sgn","Atom_name","Atom_num",...
    "sites","symmetry_operation","klist_cart","klist_l","klist_frac","kpoints_l","kpoints_name",...
    "Rnn","nn_store"...
    ];
    for i = 1:numel(CopyItem)
        TBkitobj.(CopyItem(i)) = TBkitobj_in.(CopyItem(i));
    end
end

