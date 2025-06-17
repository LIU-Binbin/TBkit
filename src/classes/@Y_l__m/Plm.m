% The associated Legendre polynomials Plm can be defined recursively as:
function PlmExpr = Plm(l,m,seed,options)
% \begin{equation}
% \begin{aligned}
% P_{0}^{0}(\cos \theta) &=1 \\
% P_{m}^{m}(\cos \theta) &=(-2l+1) \sin \theta P_{m-1}^{m-1}(\cos \theta) \\
% P_{l=m+1}^{m}(\cos \theta) &=(2 l-1) \cos \theta P_{m}^{m}(\cos \theta) \\
% P_{l}^{m}(\cos \theta) &=\frac{2 l-1}{l-m} \cos \theta P_{l-1}^{m}(\cos \theta)+\frac{1-l-m}{l-m} P_{l-2}^{m}(\cos \theta)
% \end{aligned}
% \end{equation}
%
arguments
    l double{mustBeNonnegative,mustBeInteger};
    m double{mustBeInteger}
    seed sym{} = (sym('theta','real'));
    options.ClosedForm = false;
    options.triangle = true;
end
optionsCell = namedargs2cell(options);
if abs(m)>l
    ME = MException('vasplib:Y_l__m:WrongInput', ...
        'Variable m:%d out of the range of [-l,l](l:%d)',m,l);
    throw(ME);
end
if strcmp(string(seed),"theta") && ~options.triangle
    x = sym('x','real');
    y = (1-x^2)^0.5;
elseif options.triangle
    x  = cos(seed);
    y  = sin(seed);
else
    x = seed;
    y = (1-x^2)^0.5;
end
if options.ClosedForm
    %
    % \begin{equation}
    % P_{l}^{m}(x)=(-1)^{m} \cdot 2^{l} \cdot\left(1-x^{2}\right)^{m / 2} \cdot \sum_{k=m}^{l} \frac{k !}{(k-m) !} \cdot x^{k-m} \cdot\left(\begin{array}{l}
    % l \\
    % k
    % \end{array}\right)\left(\begin{array}{l}
    % \frac{l+k-1}{2} \\
    % l
    % \end{array}\right)
    % \end{equation}
    if  options.triangle
        PlmExpr_pre = (-1)^m * 2^l * (y)^m;
    else
        PlmExpr_pre = (-1)^m * 2^l * ((1-x^2)^(1/2))^m;
    end
    % binomialfunc = @Y_l__m.binomial
    binomialfunc = @nchoosek;
    k = m;
    PlmExpr_tail = factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
    for k = m+1:l
        PlmExpr_tail =PlmExpr_tail + factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
    end
    PlmExpr = simplify(PlmExpr_pre*simplify(PlmExpr_tail));
    return;
end
% P_{\ell}^{-m}=(-1)^{m} \frac{(\ell-m) !}{(\ell+m) !} P_{\ell}^{m}
if m < 0
    PlmExpr = (-1)^(-m) * factorial(l+m)/factorial(l-m) * Y_l__m.Plm(-m,l,seed,optionsCell{:});
    return;
end
if l==0 && m ==0
    PlmExpr = 1;
    return;
elseif l == m
    % first rule
    %PlmExpr = -(2*l-1)*y*Y_l__m.Plm(l-1,m-1,seed,optionsCell{:});
    PlmExpr = (-1)^l*Y_l__m.DFactorial(2*l-1)*y^l;
    return;
elseif m == l-1
    % second rule
    PlmExpr = (2*m+1)*x*Y_l__m.Plm(m,m,seed,optionsCell{:});
else
    PlmExpr = (2*l-1)/(l-m) *x * Y_l__m.Plm(l-1,m,seed,optionsCell{:})...
        + (1-l-m)/(l-m) * Y_l__m.Plm(l-2,m,seed,optionsCell{:});
end
end
% fully normalized associated Legendre polynomials Nlm
