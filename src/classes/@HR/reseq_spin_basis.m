function H_hr = reseq_spin_basis(H_hr, old2new)
arguments
    H_hr HR
    old2new {mustBeMember(old2new,{'uudd2udud','udud2uudd'})}
end
WAN_NUM = H_hr.WAN_NUM;
if old2new == "uudd2udud"
    H_hr = H_hr.reseq(reshape([1:(WAN_NUM/2); (1:(WAN_NUM/2))+WAN_NUM/2],1, WAN_NUM));

elseif old2new == "udud2uudd"
    H_hr = H_hr.reseq([1:2:(WAN_NUM-1),2:2:(WAN_NUM)]);

end

end