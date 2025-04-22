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
                tempOrbSites = repmat(sites(i),[1,2]);
                tempOrbSites(1).orb_sym = "s";tempOrbSites(1).orb = "s";
                tempOrbSites(2).orb_sym = "z";tempOrbSites(2).orb = "p";
            case 'pxpy'
                tempOrbSites = repmat(sites(i),[1,2]);
                tempOrbSites(1).orb_sym = "x";tempOrbSites(1).orb = "p";
                tempOrbSites(2).orb_sym = "y";tempOrbSites(2).orb = "p";
            case 'p'
                tempOrbSites = repmat(sites(i),[1,3]);
                tempOrbSites(1).orb_sym = "z";tempOrbSites(1).orb = "p";
                tempOrbSites(2).orb_sym = "x";tempOrbSites(2).orb = "p";
                tempOrbSites(3).orb_sym = "y";tempOrbSites(3).orb = "p";
            case 'sd'
                tempOrbSites = repmat(sites(i),[1,6]);
                tempOrbSites(1).orb_sym = "s";tempOrbSites(1).orb = "s";
                tempOrbSites(2).orb_sym = "z^2";tempOrbSites(2).orb = "d";
                tempOrbSites(3).orb_sym = "xz";tempOrbSites(3).orb = "d";
                tempOrbSites(4).orb_sym = "yz";tempOrbSites(4).orb = "d";
                tempOrbSites(5).orb_sym = "x^2-y^2";tempOrbSites(5).orb = "d";
                tempOrbSites(6).orb_sym = "xy";tempOrbSites(6).orb = "d";
            case 'spd'
                tempOrbSites = repmat(sites(i),[1,9]);
                tempOrbSites(1).orb_sym = "s";tempOrbSites(1).orb = "s";
                tempOrbSites(2).orb_sym = "z";tempOrbSites(2).orb = "p";
                tempOrbSites(3).orb_sym = "x";tempOrbSites(3).orb = "p";
                tempOrbSites(4).orb_sym = "y";tempOrbSites(4).orb = "p";
                tempOrbSites(5).orb_sym = "z^2";tempOrbSites(5).orb = "d";
                tempOrbSites(6).orb_sym = "xz";tempOrbSites(6).orb = "d";
                tempOrbSites(7).orb_sym = "yz";tempOrbSites(7).orb = "d";
                tempOrbSites(8).orb_sym = "x^2-y^2";tempOrbSites(8).orb = "d";
                tempOrbSites(9).orb_sym = "xy";tempOrbSites(9).orb = "d";
            case 'pd'
                tempOrbSites = repmat(sites(i),[1,8]);
                tempOrbSites(1).orb_sym = "z";tempOrbSites(1).orb = "p";
                tempOrbSites(2).orb_sym = "x";tempOrbSites(2).orb = "p";
                tempOrbSites(3).orb_sym = "y";tempOrbSites(3).orb = "p";
                tempOrbSites(4).orb_sym = "z^2";tempOrbSites(4).orb = "d";
                tempOrbSites(5).orb_sym = "xz";tempOrbSites(5).orb = "d";
                tempOrbSites(6).orb_sym = "yz";tempOrbSites(6).orb = "d";
                tempOrbSites(7).orb_sym = "x^2-y^2";tempOrbSites(7).orb = "d";
                tempOrbSites(8).orb_sym = "xy";tempOrbSites(8).orb = "d";
            case 'd'
                tempOrbSites = repmat(sites(i),[1,5]);
                tempOrbSites(1).orb_sym = "z^2";tempOrbSites(1).orb = "d";
                tempOrbSites(2).orb_sym = "xz";tempOrbSites(2).orb = "d";
                tempOrbSites(3).orb_sym = "yz";tempOrbSites(3).orb = "d";
                tempOrbSites(4).orb_sym = "x^2-y^2";tempOrbSites(4).orb = "d";
                tempOrbSites(5).orb_sym = "xy";tempOrbSites(5).orb = "d";
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