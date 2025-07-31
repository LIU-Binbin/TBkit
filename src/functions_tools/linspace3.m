function vlist = linspace3(v1, v2, n)
dim = length(v1);
vlist = zeros(n, dim);
for i = 1:dim
    vlist(1:n,i) = linspace(v1(i),v2(i),n);
end
end