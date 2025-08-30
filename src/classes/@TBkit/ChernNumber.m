function Chern = ChernNumber(Ham, tensor_index, optionsK, options)
arguments
    Ham TBkit;
    tensor_index (1,2) double = [1 2] % Omega_xy

    optionsK.kstart(1,3) double = [0 0 0]
    optionsK.kdir1 (1,3) double = [1 0 0]
    optionsK.kdir2 (1,3) double = [0 1 0]
    optionsK.Nk1 double = 50
    optionsK.Nk2 double = 50

    options.BAND_index = [];
end
optionsKcell = namedargs2cell(optionsK);
optionsCell = namedargs2cell(options);

BCCAR = BerryCurvature_2D(Ham, tensor_index, optionsKcell{:},optionsCell{:});

Ksquare = norm(cross(optionsK.kdir1*Ham.Gk, optionsK.kdir2*Ham.Gk));
Chern = - sum(BCCAR,'all')/(2*pi) /optionsK.Nk1/optionsK.Nk2 * Ksquare;
Chern = roundn(Chern, round(log10(1e-3)));
end