function str_out = quadk2K(str_in)
str_out = strrep(str_in,'x^2',"X");
str_out = strrep(str_out,'y^2',"Y");
str_out = strrep(str_out,'z^2',"Z");
str_out = strrep(str_out,'x^3',"X*x");
str_out = strrep(str_out,'y^3',"Y*y");
str_out = strrep(str_out,'z^3',"Z*z");
str_out = strrep(str_out,'x^4',"X^2");
str_out = strrep(str_out,'y^4',"Y^2");
str_out = strrep(str_out,'z^4',"Z^2");
str_out = strrep(str_out,'x^5',"X^2*x");
str_out = strrep(str_out,'y^5',"Y^2*y");
str_out = strrep(str_out,'z^5',"Z^2*z");
end
