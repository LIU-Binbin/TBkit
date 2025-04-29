function R = Sph2Cart(S)
            % R =S;
            r= S(:,1);
            theta = S(:,2);
            phi = S(:,3);
            z = r .*cos(theta);
            rcoselev = r .* sin(theta);
            x = rcoselev .* cos(phi);
            y = rcoselev .* sin(phi);
            R = [x,y,z];
        end