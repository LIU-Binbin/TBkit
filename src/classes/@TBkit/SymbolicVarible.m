function Symble = SymbolicVarible(seeds,superscript,subscript,level)
Superscript = "";
Subscript = "";
for i = 1:length(superscript)
    if superscript(i) < 0
        Superscript = Superscript+"__"+num2str(round(abs(superscript(i))))+"_bar";
    else
        Superscript = Superscript+"__"+num2str(round(abs(superscript(i))));
    end
end
for i = 1:length(subscript)
    if subscript(i) < 0
        Subscript = Subscript+"_"+num2str(round(abs(subscript(i))))+"_bar";
    else
        Subscript = Subscript+"_"+num2str(round(abs(subscript(i))));
    end
end
if nargin == 4
    Level = "_"+num2str(level)+"_ubar";
    Symble = sym(seeds+Subscript+Superscript+Level,'real');
else
    Symble = sym(seeds+Subscript+Superscript,'real');
end
end