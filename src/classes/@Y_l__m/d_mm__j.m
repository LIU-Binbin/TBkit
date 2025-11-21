function SymExpr = d_mm__j(j,m1,m2,theta)
% d_mm__j  ---  Jacobi polynomial formula for Wigner small-d
% Implements the Rose/Edmonds region formula:
%
%   d^j_{m2,m1}(theta)
%   = (-1)^lambda * sqrt(C1/C2)
%     * (sin(theta/2))^a (cos(theta/2))^b
%     * P_k^{(a,b)}(cos theta)
%
% Regions determined by comparing:
%   R1 = j+m2
%   R2 = j-m2
%   R3 = j+m1
%   R4 = j-m1
%
% The correct k = min(R1,R2,R3,R4)
%
% NOTE: This version preserves your structure,
%       but fixes the missing symmetry factors and ambiguous region ties.

arguments
    j double {mustBePositive}
    m1 double
    m2 double
    theta sym = sym('theta')
end

% -------------------------
% enforce valid domain
% -------------------------
if abs(m1) > j || abs(m2) > j
    SymExpr = sym(0);
    return;
end

% -------------------------
% symmetry: Rose  (must include (-1)^{m2-m1})
% -------------------------
if m1 < 0 && m2 < 0
    % Wigner d symmetry:
    % d^j_{-m2,-m1}(θ) = (-1)^(m1-m2) d^j_{m2,m1}(θ)
    SymExpr = (-1)^(m1-m2) * Y_l__m.d_mm__j(j,-m1,-m2,theta);
    return;
end

% -------------------------
% region numbers (Rose)
% -------------------------
R1 = j+m2;
R2 = j-m2;
R3 = j+m1;
R4 = j-m1;

k = min([R1,R2,R3,R4]);

% -------------------------
% determine the correct region
% Rose's rule:
%   If k == R1 → region 1
%   If k == R2 → region 2
%   If k == R3 → region 3
%   If k == R4 → region 4
% (multiple matches need tie-breaking based on m1,m2 signs)
% -------------------------

if k == R1 && ~(k==R2 || k==R3 || k==R4)
    region = 1;
elseif k == R2 && ~(k==R3 || k==R4)
    region = 2;
elseif k == R3 && ~(k==R4)
    region = 3;
else
    region = 4;
end

% -------------------------
% assign a,b,lambda based on region
% -------------------------
switch region
    case 1    % R1 = j+m2 minimal
        a = m1 - m2;
        lambda = m1 - m2;

    case 2    % R2 = j-m2 minimal
        a = m2 - m1;
        lambda = 0;

    case 3    % R3 = j+m1 minimal
        a = m2 - m1;
        lambda = 0;

    case 4    % R4 = j-m1 minimal
        a = m1 - m2;
        lambda = m1 - m2;
end

% -------------------------
% b from Rose's relation
% -------------------------
b = 2*j - 2*k - a;

% -------------------------
% Jacobi Polynomial factor
% -------------------------
SymExpr = (-1)^lambda * ...
    sqrt(nchoosek(2*j-k, k+a) / nchoosek(k+b, b)) * ...
    sin(theta/2)^a * cos(theta/2)^b * ...
    jacobiP(k, a, b, cos(theta));

end
