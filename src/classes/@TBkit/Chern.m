function Chern_number = Chern(TBkitobj,options,options_Chern)
%CHERN Calculate Chern number for 2D system
%
%   Syntax:
%       Chern_number = Chern(TBkitobj,options,options_Chern)
%
%   Description:
%       Computes Chern number by integrating Berry curvature over BZ.
%
%   Inputs:
%       TBkitobj     - Tight-binding model object
%       options      - k-mesh options:
%                      knum1/2   - Grid dimensions
%                      kstart    - Starting k-point
%                      kdir1/2   - Grid directions
%                      cartesian - Cartesian coordinates flag
%                      dir_seq   - Direction sequence
%                      dir_start - Starting direction
%       options_Chern - Calculation options:
%                      Accuracy  - Numerical tolerance
%
%   Output:
%       Chern_number - Integer Chern number
arguments
    TBkitobj ;
    options.BAND_index = [];
    options.knum1   = 51;
    options.knum2   = 51;
    options.kstart  = [-0.5,-0.5,0];
    options.kdir1   = [1,0,0];
    options.kdir2   = [0,1,0];
    options.cartesian = false;
    options.dir_seq = [1,2,3];
    options.dir_start = 'k_z';
    options.fig = false;
    options.plot = false;
    options.Oper = [];
    options.subband = [];
    options.ProjectionMethod = 'sign';
    options.ProjectionStruct = struct('field','imag');
    options_Chern.Accuracy = 1e-6;
end
optionsCell = namedargs2cell(options);
[BCCAR,~,~,~] = BC_2D(TBkitobj,optionsCell{:});
Chern_number = roundn(sum(BCCAR,'all')/(2*pi),round(log10(options_Chern.Accuracy)));
end