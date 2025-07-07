classdef constants
    %% physical constants
    properties (Constant)
        h_eV_s = 4.1357e-15; % eV*s
        hbar_eV_s = 6.5821e-16; % eV*s
        charge_C = 1.6021766e-19; % Coulomb
        muB_eV_T = 5.7884e-5; % eV/Tesla
        kB_eV_K = 8.61733e-5; % eV/Kelvin
    end
    %% matlab_default_colors
    properties (Constant)
        blue_dark = [0 0.4470 0.7410];
        orange = [0.8500 0.3250 0.0980];
        yellow = [0.9290 0.6940 0.1250];
        purple = [0.4940 0.1840 0.5560];
        green = [0.4660 0.6740 0.1880];
        blue_light = [0.3010 0.7450 0.9330];
        red = [0.6350 0.0780 0.1840];
    end

    methods
        function obj = constants()
        end
    end
end