function TBkitobj = TBkit_init(TBkitobj)
    if isempty(TBkitobj.Rm)
        TBkitobj.Rm = eye(TBkitobj.Dim);
    end
    if isempty(TBkitobj.VarsSeqLcart)
        syms k_x k_y k_z k_w real;
        TBkitobj.VarsSeqLcart = [k_x k_y k_z k_w];
    end
    if isempty(TBkitobj.VarsSeqLfrac)
        syms k_1 k_2 k_3 k_4 real;
        TBkitobj.VarsSeqLfrac = [k_1 k_2 k_3 k_4];
    end
end
