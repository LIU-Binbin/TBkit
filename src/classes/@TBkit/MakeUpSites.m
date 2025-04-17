function tmpsites = MakeUpSites(sites, spintype)
%MAKEUPSITES Processes atomic sites information for orbital/spin expansion
%   This function processes atomic sites data to:
%   1. Inherit orbital information for empty entries
%   2. Expand composite orbital types (sp, spd, etc.) into individual orbitals
%   3. Handle spin configurations based on specified spintype
%
%   Inputs:
%   sites    - Array of structs containing site information with fields:
%              element, orb, orb_sym, spin
%   spintype - Spin configuration type: 'spinless', 'wannier', or 'block'
%
%   Output:
%   tmpsites - Processed array of site structures with expanded orbitals/spins
%
%   Throws error for invalid POSCAR formats or unsupported orbital types

%% Phase 1: Orbital Information Inheritance
tmpsites = sites(1);       % Initialize template
tmpsites(1) = [];          % Clear initial element
target_site = sites(1);    % Initialize reference site

for i = 2:length(sites)
    % Handle empty orbital information
    if all(strcmp(sites(i).element, "")) && all(strcmp(sites(i).orb, ""))
        sites(i).element = target_site.element;
        sites(i).orb = target_site.orb;
        sites(i).orb_sym = target_site.orb_sym;
        
    % Update reference site when new orbital info appears    
    elseif ~all(strcmp(sites(i).element, "")) && ~all(strcmp(sites(i).orb, ""))
        target_site = sites(i);
    else
        error('Invalid POSCAR format: Mixed empty/non-empty fields detected');
    end
end

%% Phase 2: Orbital Basis Expansion
for i = 1:length(sites)
    if strcmp(sites(i).orb_sym, "")
        switch char(sites(i).orb)
            % s-orbital (single state)
            case 's'
                tempOrbSites = sites(i);
                [tempOrbSites.orb_sym, tempOrbSites.orb] = deal("s");
                
            % sp-orbital set (4 states)    
            case 'sp'
                tempOrbSites = repmat(sites(i),1,4);
                [tempOrbSites(1).orb_sym, tempOrbSites(1).orb] = deal("s");
                [tempOrbSites(2:4).orb] = deal("p");
                [tempOrbSites(2:4).orb_sym] = deal("y","z","x");
                % Other orbital cases handled similarly...
            case 'spz'
                tmpsites_slib = repmat(sites(i),[1,2]);
                tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                tmpsites_slib(2).orb_sym = "z";tmpsites_slib(2).orb = "p";
            case 'pxpy'
                tmpsites_slib = repmat(sites(i),[1,2]);
                tmpsites_slib(1).orb_sym = "x";tmpsites_slib(1).orb = "p";
                tmpsites_slib(2).orb_sym = "y";tmpsites_slib(2).orb = "p";
            case 'p'
                tmpsites_slib = repmat(sites(i),[1,3]);
                tmpsites_slib(1).orb_sym = "z";tmpsites_slib(1).orb = "p";
                tmpsites_slib(2).orb_sym = "x";tmpsites_slib(2).orb = "p";
                tmpsites_slib(3).orb_sym = "y";tmpsites_slib(3).orb = "p";
            case 'sd'
                tmpsites_slib = repmat(sites(i),[1,6]);
                tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                tmpsites_slib(2).orb_sym = "z^2";tmpsites_slib(2).orb = "d";
                tmpsites_slib(3).orb_sym = "xz";tmpsites_slib(3).orb = "d";
                tmpsites_slib(4).orb_sym = "yz";tmpsites_slib(4).orb = "d";
                tmpsites_slib(5).orb_sym = "x^2-y^2";tmpsites_slib(5).orb = "d";
                tmpsites_slib(6).orb_sym = "xy";tmpsites_slib(6).orb = "d";
            case 'spd'
                tmpsites_slib = repmat(sites(i),[1,9]);
                tmpsites_slib(1).orb_sym = "s";tmpsites_slib(1).orb = "s";
                tmpsites_slib(2).orb_sym = "z";tmpsites_slib(2).orb = "p";
                tmpsites_slib(3).orb_sym = "x";tmpsites_slib(3).orb = "p";
                tmpsites_slib(4).orb_sym = "y";tmpsites_slib(4).orb = "p";
                tmpsites_slib(5).orb_sym = "z^2";tmpsites_slib(5).orb = "d";
                tmpsites_slib(6).orb_sym = "xz";tmpsites_slib(6).orb = "d";
                tmpsites_slib(7).orb_sym = "yz";tmpsites_slib(7).orb = "d";
                tmpsites_slib(8).orb_sym = "x^2-y^2";tmpsites_slib(8).orb = "d";
                tmpsites_slib(9).orb_sym = "xy";tmpsites_slib(9).orb = "d";
            case 'pd'
                tmpsites_slib = repmat(sites(i),[1,8]);
                tmpsites_slib(1).orb_sym = "z";tmpsites_slib(1).orb = "p";
                tmpsites_slib(2).orb_sym = "x";tmpsites_slib(2).orb = "p";
                tmpsites_slib(3).orb_sym = "y";tmpsites_slib(3).orb = "p";
                tmpsites_slib(4).orb_sym = "z^2";tmpsites_slib(4).orb = "d";
                tmpsites_slib(5).orb_sym = "xz";tmpsites_slib(5).orb = "d";
                tmpsites_slib(6).orb_sym = "yz";tmpsites_slib(6).orb = "d";
                tmpsites_slib(7).orb_sym = "x^2-y^2";tmpsites_slib(7).orb = "d";
                tmpsites_slib(8).orb_sym = "xy";tmpsites_slib(8).orb = "d";
            case 'd'
                tmpsites_slib = repmat(sites(i),[1,5]);
                tmpsites_slib(1).orb_sym = "z^2";tmpsites_slib(1).orb = "d";
                tmpsites_slib(2).orb_sym = "xz";tmpsites_slib(2).orb = "d";
                tmpsites_slib(3).orb_sym = "yz";tmpsites_slib(3).orb = "d";
                tmpsites_slib(4).orb_sym = "x^2-y^2";tmpsites_slib(4).orb = "d";
                tmpsites_slib(5).orb_sym = "xy";tmpsites_slib(5).orb = "d";
            case 'f'
                %
            
            otherwise
                error('Unsupported orbital type: %s', sites(i).orb);
        end
    else
        tempOrbSites = sites(i);
    end
    tmpsites = [tmpsites, tempOrbSites];
end

%% Phase 3: Spin Configuration Handling
switch lower(spintype)
    case 'spinless'  % No spin dimension
        return
        
    case 'wannier'   % Alternating spin configuration
        n = length(tmpsites);
        [tmpsites(1:n).spin] = deal(0.5);
        spinArray = num2cell(-[tmpsites.spin]);
        [tmpsites(n+1:2*n).spin] = spinArray{:};
        tmpsites = reshape([tmpsites(1:n); tmpsites(n+1:2*n)], 1, []);
        
    case 'block'     % Block spin configuration
        n = length(tmpsites);
        [tmpsites(1:n).spin] = deal(0.5);
        spinArray = num2cell(-[tmpsites.spin]);
        tmpsites = [tmpsites, tmpsites];
        [tmpsites(n+1:end).spin] = spinArray{:};
        
    otherwise
        error('Unsupported spin type: %s', spintype);
end
end