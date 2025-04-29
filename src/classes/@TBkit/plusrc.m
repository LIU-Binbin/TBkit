function rc_plus = plusrc(rc)
%PLUSRC Adjust real coordinate to be within [0,1) range
%   RC_PLUS = PLUSRC(RC) adjusts input coordinate RC to be within the
%   fundamental domain [0,1) by:
%     - Adding 1 if RC < 0
%     - Taking modulo 1 if RC >= 1
%     - Returning as-is if 0 <= RC < 1
%
%   Input:
%       rc - Input coordinate value
%
%   Output:
%       rc_plus - Adjusted coordinate in [0,1) range
if rc<0
    %disp('bugtest');
    rc_plus = rc +1;

elseif rc>=1
    rc_plus =mod(rc,1);
else
    rc_plus =rc;
end
end