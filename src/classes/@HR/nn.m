function H_hr = nn(H_hr,search_range,Accuracy,Rlength_cut,options)
arguments
H_hr HR;
search_range double = [0 0 0];
Accuracy double = 1e-4;
Rlength_cut double = 5;
options.MAX_NN_LENGTH =  10000000;
options.onsite = false;
end
optionsCell = namedargs2cell(options);
H_hr = nn@TBkit(H_hr,search_range,Accuracy,Rlength_cut,optionsCell{:});
end
