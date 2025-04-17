        function rc_plus = plusrc(rc)
            if rc<0
                %disp('bugtest');
                rc_plus = rc +1;

            elseif rc>=1
                rc_plus =mod(rc,1);
            else
                rc_plus =rc;
            end
        end