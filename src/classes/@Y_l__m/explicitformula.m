function expression = explicitformula(YlmObj,options)
arguments
YlmObj Y_l__m;
options.vpa = true;
options.explicit = true;
options.cart = false
options.seed =["theta","phi"];
end
optionsCell = namedargs2cell(options);
if YlmObj.m < 0
YlmObj.m = -YlmObj.m;
expression = conj(explicitformula(YlmObj,optionsCell{:}))*(-1)^YlmObj.m;
return;
end
if ~options.cart
theta = sym(options.seed(1),'real');
phi = sym(options.seed(2),'real');
else
syms x y z r real;
end
if options.cart
coe_pi = sqrt(1/sym(pi));
coe_m = (-1)^YlmObj.m;
coe_r = 1/r^YlmObj.l;
XIYfunc = (x+1i*y)^YlmObj.m;
switch  YlmObj.l
case 0
coe_1 = 1/2;
coe_front_pi = 1;
Zfunc = 1;
case 1
coe_1 = 1/2;
switch YlmObj.m
case 0
coe_front_pi = sqrt(3);
Zfunc = z;
case 1
coe_front_pi = sqrt(3/2);
Zfunc = 1;
end
case 2
coe_1 = 1/4;
switch YlmObj.m
case 0
coe_front_pi = sqrt(5);
Zfunc = 3*z^2-r^2;
case 1
coe_front_pi = sqrt(30);
Zfunc = z;
case 2
coe_front_pi = sqrt(15/2);
Zfunc = 1;
end
case 3
coe_1 = 1/4;
switch YlmObj.m
case 0
coe_front_pi = sqrt(7);
Zfunc = (5*z^2-3*r^2)*z;
case 1
coe_front_pi = sqrt(21/4);
Zfunc = (5*z^2-r^2);
case 2
coe_front_pi = sqrt(105/2);
Zfunc = z;
case 3
coe_front_pi = sqrt(35/4);
Zfunc = 1;
end
case 4
coe_1 = 3/16;
switch YlmObj.m
case 0
coe_front_pi = sqrt(1);
Zfunc = (35*z^4-30*z^2*r^2+3*r^2);
case 1
coe_front_pi = sqrt(40);
Zfunc = (7*z^2-3*r^2)*z;
case 2
coe_front_pi = sqrt(8);
Zfunc = (7*z^2-r^2);
case 3
coe_front_pi = sqrt(35*4);
Zfunc = z;
case 4
coe_front_pi = sqrt(35/2);
Zfunc = 1;
end
case 5
coe_1 = 1/16;
switch YlmObj.m
case 0
coe_front_pi = sqrt(11);
Tfunc = SINT^YlmObj.m*(63*COST^5-70*COST^3+15*COST);
case 1
coe_front_pi = sqrt(165/2);
Tfunc = SINT^YlmObj.m*(21*COST^4-14*COST^2+1);
case 2
coe_front_pi = sqrt(1155*2);
Tfunc = SINT^YlmObj.m*(3*COST^3-COST);
case 3
coe_front_pi = sqrt(385/4);
Tfunc = SINT^YlmObj.m*(9*COST^2-1);
case 4
coe_front_pi = sqrt(385*9/2);
Tfunc = SINT^YlmObj.m*COST;
case 5
coe_front_pi = sqrt(77*9/4);
Tfunc = SINT^YlmObj.m;
end
otherwise
expression =Y_l__m.InnerY_l__m(YlmObj.l,YlmObj.m,'triangle',false);
return;
end
expression = coe_m*coe_1*coe_front_pi*coe_pi*XIYfunc*Zfunc*coe_r;
else
EIPHI = exp(1i*YlmObj.m*phi);
COST = cos(theta);SINT = sin(theta);
coe_pi = sqrt(1/sym(pi));
coe_m = (-1)^YlmObj.m;
switch  YlmObj.l
case 0
coe_1 = 1/2;
coe_front_pi = 1;
Tfunc = 1;
case 1
coe_1 = 1/2;
switch YlmObj.m
case 0
coe_front_pi = sqrt(3);
Tfunc = COST;
case 1
coe_front_pi = sqrt(3/2);
Tfunc = SINT;
end
case 2
coe_1 = 1/4;
switch YlmObj.m
case 0
coe_front_pi = sqrt(5);
Tfunc = 3*COST^2-1;
case 1
coe_front_pi = sqrt(30);
Tfunc = SINT*COST;
case 2
coe_front_pi = sqrt(15/2);
Tfunc = SINT^2;
end
case 3
coe_1 = 1/4;
switch YlmObj.m
case 0
coe_front_pi = sqrt(7);
Tfunc = SINT^YlmObj.m*(5*COST^3-3*COST);
case 1
coe_front_pi = sqrt(21/4);
Tfunc = SINT^YlmObj.m*(5*COST^2-1);
case 2
coe_front_pi = sqrt(105/2);
Tfunc = SINT^YlmObj.m*COST;
case 3
coe_front_pi = sqrt(35/4);
Tfunc = SINT^YlmObj.m;
end
case 4
coe_1 = 3/16;
switch YlmObj.m
case 0
coe_front_pi = sqrt(1);
Tfunc = SINT^YlmObj.m*(35*COST^4-30*COST^2+3);
case 1
coe_front_pi = sqrt(40);
Tfunc = SINT^YlmObj.m*(7*COST^3-3*COST);
case 2
coe_front_pi = sqrt(8);
Tfunc = SINT^YlmObj.m*(7*COST^2-1);
case 3
coe_front_pi = sqrt(35*4);
Tfunc = SINT^YlmObj.m*COST;
case 4
coe_front_pi = sqrt(35/2);
Tfunc = SINT^YlmObj.m;
end
case 5
coe_1 = 1/16;
switch YlmObj.m
case 0
coe_front_pi = sqrt(11);
Tfunc = SINT^YlmObj.m*(63*COST^5-70*COST^3+15*COST);
case 1
coe_front_pi = sqrt(165/2);
Tfunc = SINT^YlmObj.m*(21*COST^4-14*COST^2+1);
case 2
coe_front_pi = sqrt(1155*2);
Tfunc = SINT^YlmObj.m*(3*COST^3-COST);
case 3
coe_front_pi = sqrt(385/4);
Tfunc = SINT^YlmObj.m*(9*COST^2-1);
case 4
coe_front_pi = sqrt(385*9/2);
Tfunc = SINT^YlmObj.m*COST;
case 5
coe_front_pi = sqrt(77*9/4);
Tfunc = SINT^YlmObj.m;
end
otherwise
expression =Y_l__m.InnerY_l__m(YlmObj.l,YlmObj.m);
return;
end
expression = coe_m*coe_1*coe_front_pi*coe_pi*Tfunc*EIPHI;
end
end
