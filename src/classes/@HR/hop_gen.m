function hop = hop_gen(hop_pre,nn_level)
tempstr = string(simplify(hop_pre));
tempstr = strcat(tempstr,"_",string(nn_level));
hop = sym(tempstr,'real');
end
