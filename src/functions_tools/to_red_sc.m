function To_red_sc = to_red_sc(red_vec_orig, Ns)
    % Compute To_red_sc by multiplying red_vec_orig with the inverse of Ns
    % using the matrix division operator for efficiency.
    To_red_sc = red_vec_orig / Ns;
end
