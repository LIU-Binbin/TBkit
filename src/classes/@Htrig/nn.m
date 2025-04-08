function H_htrig = nn(H_htrig,search_range,Accuracy,Rlength_cut)
if nargin <4
Rlength_cut = 5;
end
if nargin <3
Accuracy = 4;
end
if nargin <2
search_range = [0 0 0];
end
H_htrig = nn@TBkit(H_htrig,search_range,Accuracy,Rlength_cut);
H_htrig.Type = 'exp';
end
