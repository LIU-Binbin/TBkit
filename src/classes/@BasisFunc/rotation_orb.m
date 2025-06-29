function BForb = rotation_orb(BForb, Rf, tf, options)
%ROTATION_ORB  Apply rotation and translation to orbital positions.
%
%   BForb = ROTATION_ORB(BForb, Rf, tf, options) rotates and translates the orbital
%   position vector `BForb` by applying a rotation matrix `Rf` and a translation
%   vector `tf`, with optional behavior specified via `options`.
%
%   Inputs:
%       BForb              - 1x3 orbital position vector (or compatible array).
%       Rf                 - 3x3 rotation matrix applied to orbital positions 
%                            (default: identity).
%       tf                 - 3x1 translation vector in fractional coordinates 
%                            (default: [0 0 0]).
%       options            - Structure with optional arguments (currently mostly ignored):
%           sym            - Whether to use symbolic simplification (default: true).
%           conjugate      - Whether to apply complex conjugation (default: false).
%           antisymmetry   - Whether antisymmetry affects result (default: false).
%           forgetcoe      - Whether to ignore coefficients (default: false).
%           fast           - Toggle for fast implementation (default: true).
%           hybird         - Not implemented hybrid mode (default: false).
%           spincoupled    - Whether spin is coupled in transformation (default: false).
%           orbcoupled     - Whether orbital is coupled (default: false).
%           raw            - Toggle for raw representation (default: true).
%           vpalevel       - VPA precision level if symbolic (default: 6).
%           center         - Center for rotation (default: [0,0,0]).
%
%   Output:
%       BForb              - Rotated and translated orbital position (wrapped to [0,1)^3).
%
%   Behavior:
%       1. Applies affine transformation:
%          BForb_new = (BForb - center) * Rf.' + center + tf
%       2. Each coordinate is then wrapped back into unit cell via mod 1.
%
%   Example:
%       BForb = [0.25, 0.5, 0.75];
%       Rf = rotation_matrix_z(pi);
%       tf = [0; 0; 0];
%       BForb_new = rotation_orb(BForb, Rf, tf);
%
%   Notes:
%       - If `BForb` is empty, function returns immediately.
%       - All entries are assumed to be in **fractional coordinates**.
%       - `options` fields are parsed for compatibility but not used actively.
%
%   See also: rotate, mod

    arguments
        BForb
        Rf {mustBeSize(Rf,[3 3])} = diag([1 1 1]);  % Rotation matrix
        tf {mustBeSize(tf,[3 1;1 3])} = [0 0 0];     % Translation vector
        options.sym = true;
        options.conjugate = false;
        options.antisymmetry = false;
        options.forgetcoe = false;
        options.fast = true;
        options.hybird = false;
        options.spincoupled = false;
        options.orbcoupled = false;
        options.raw = true;
        options.vpalevel = 5;
        options.center = [0, 0, 0];
    end

    if isempty(BForb)
        return;  
    end

    BForb = roundn( double((BForb - options.center) * Rf.' + options.center + tf),-options.vpalevel);

    for i = 1:3
        BForb(i) = mod(BForb(i), 1);
    end
end

