function BForb = rotation_orb(BForb,Rf,tf,options)
arguments
BForb
Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);
tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);
options.sym = true;
options.conjugate = false;
options.antisymmetry = false;
options.forgetcoe = false;
options.fast = true;
options.hybird = false;
options.spincoupled = false;
options.orbcoupled = false;
options.raw = true;
options.vpalevel = 6;
options.center = [0,0,0];
end
if isempty(BForb)
return;
end
BForb = (BForb-options.center)*Rf.'+options.center + tf;
for i = 1:3
BForb(i) =mod(BForb(i),1);
end
end
