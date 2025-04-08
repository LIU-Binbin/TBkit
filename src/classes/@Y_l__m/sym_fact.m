function sym_fact = sym_fact(n)
if n == 0
sym_fact = sym(1);
end
if n > 0
sym_fact = prod(sym(1):sym(n));
end
if n < 0
sym_fact = sym(inf);
end
end
