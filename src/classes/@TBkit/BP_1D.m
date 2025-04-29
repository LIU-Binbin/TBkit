function [BP,WAVECAR_loop] = BP_1D(TBkitobj,klist,options,optionsWAVE)
%BP_1D Calculate Berry phase in 1D systems
%
%   Syntax:
%       [BP,WAVECAR_loop] = BP_1D(TBkitobj,klist,options,optionsWAVE)
%
%   Description:
%       Computes Berry phase (polarization) for a 1D band structure
%       using Wilson loop method.
%
%   Inputs:
%       TBkitobj   - Tight-binding model object
%       klist      - Array of k-points along the path
%       options    - Calculation options:
%                    BAND_index - Band indices to include
%                    ax         - Axis handle for plotting
%                    plot       - Visualization flag
%       optionsWAVE - Wavefunction options:
%                    LWAVE      - Save wavefunctions flag
%
%   Outputs:
%       BP          - Berry phase (polarization)
%       WAVECAR_loop - Wavefunctions along the loop
arguments
    TBkitobj
    klist
    options.BAND_index = [];
    options.ax = handle([]);
    options.plot = false;
    optionsWAVE.LWAVE = false;
end
optionsCell = namedargs2cell(options);
WAVECAR_loop = Topo1DpreWAVECAR(TBkitobj,klist,optionsCell{:});
if optionsWAVE.LWAVE

end
BP = sum(TBkit.wancenter_1D(WAVECAR_loop));
end