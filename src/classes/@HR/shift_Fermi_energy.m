function H_hr = shift_Fermi_energy(H_hr, Efermi)
if strcmp(H_hr.Type,'sparse')
    error("not implemented yet")
elseif strcmp(H_hr.Type,'list')
    error("not implemented yet")
elseif strcmp(H_hr.Type,'mat')
    if H_hr.num
        H_hr.HnumL(:,:,H_hr.Line_000) = H_hr.HnumL(:,:,H_hr.Line_000) - double(Efermi*eye(H_hr.WAN_NUM));
    end
    if H_hr.coe
        H_hr.HcoeL(:,:,H_hr.Line_000) = H_hr.HcoeL(:,:,H_hr.Line_000) - sym(   Efermi*eye(H_hr.WAN_NUM));
    end
end
end