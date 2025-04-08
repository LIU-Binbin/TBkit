classdef Basis < handle
    % Basis class to describe an orbital and perform rotations
    % The Basis function is a linear combination of:
    % Ylm, Yl^m, X, LS, j, mz, etc.

    properties
        BasisL {mustBeA(BasisL, ["BasisFunc", "cell", "double", "string", "char", "sym", "Spin", "Y_lm", "Y_l__m", "jm", "lm"])}
        BcoeL = sym([]);
        BnumL = [];
        orbL = [];
    end

    properties(Dependent)
        Basis_num;
    end

    properties(Hidden)
        Rm = [1 0 0; 0 1 0; 0 0 1];
        num = false;
        coe = true;
    end

    methods % Constructor
        function Basisobj = Basis(BasisL, BnumL, BcoeL)
            arguments
                BasisL {mustBeA(BasisL, ["BasisFunc", "cell", "double", "string", "char", "sym", "Spin", "Y_lm", "Y_l__m", "jm", "lm"])}
                BnumL double = []
                BcoeL sym = sym([])
            end
            % Handle different input types for BasisL
            switch class(BasisL)
                case 'double'
                    if isscalar(BasisL)
                        Basisobj.BasisL = BasisFunc(BasisL); % Example for BasisFunc
                    else
                        % Handle quantumL case here (not implemented in your code)
                    end
                case {'string', 'char'}
                    Basisobj = Basis.BasisDatabase(BasisL); % Lookup in database
                case 'sym'
                    Basisobj.BasisL = BasisFunc(BasisL); % Example for symbolic Basis
                case {'Spin', 'jm', 'lm'}
                    Basisobj.BasisL = BasisL;
                case {'Y_lm', 'Y_l__m'}
                    Basisobj.BasisL = BasisL;
            end
            % Default handling for BnumL and BcoeL if not provided
            Basis_num = size(Basisobj.BasisL, 1);
            if isempty(BnumL)
                Basisobj.BnumL = ones([1, Basis_num]);
            else
                Basisobj.BnumL = BnumL;
            end
            if isempty(BcoeL)
                Basisobj.BcoeL = sym(ones([1, Basis_num]));
            else
                Basisobj.BcoeL = BcoeL;
            end
        end
    end

    methods(Static) % Static Methods for Database
        function Basisobj = BasisDatabase(BasisName)
            switch BasisName
                case 'BHZ'
                    Basisobj = Basis([Spin(1/2, 1/2); Spin(3/2, 3/2, 'parity', -1); Spin(1/2, -1/2); Spin(3/2, -3/2, 'parity', -1)]);
                case 'Graphene'
                    % Add graphene case here
                case 'WSM'
                    % Add WSM case here
                otherwise
                    error('Unknown basis name: %s', BasisName);
            end
        end
    end

    methods % Getters
        function Basis_num = get.Basis_num(Basisobj)
            Basis_num = length(Basisobj.BasisL);
        end
    end

    methods % Display Method
        function disp(Basisobj)
            % Display Basis Function
            BasisFunction = Basisobj.BasisL;
            builtin('disp', Basisobj); % Call the built-in disp function
            
            % Prepare coefficient string for display
            CoeStr = arrayfun(@(x) Basisobj.getCoefficientString(x), Basisobj.BcoeL, 'UniformOutput', false);
            
            % Display Basis Function components
            switch class(BasisFunction)
                case 'BasisFunc'
                    for i = 1:size(BasisFunction, 1)
                        fprintf(CoeStr{i} + "{ " + string(BasisFunction(i).BFuncL{1}));
                        for j = 2:length(BasisFunction(i).BFuncL)
                            fprintf(' + ' + string(BasisFunction(i).BFuncL{j}));
                        end
                        fprintf(' }\n');
                    end
                case 'Spin'
                    for i = 1:size(BasisFunction, 1)
                        fprintf(CoeStr{i} + "{ " + string(BasisFunction(i, :)) + ' }\n');
                    end
            end
        end
    end

    methods % Helper Function for Coefficients
        function CoeStr = getCoefficientString(~, coeff)
            if coeff == 1
                CoeStr = "";
            elseif coeff == -1
                CoeStr = "-";
            else
                CoeStr = string(coeff);
            end
        end
    end

    methods % Rotation & Transformation
        function Umat = U(Basisobj, rotm, t, rightorleft)
            arguments
                Basisobj Basis;
                rotm;
                t = [0, 0, 0];
                rightorleft = 'right';
            end
            if isequal(size(rotm), [3, 3])
                [n, theta] = Oper.Rotation2nTheta(rotm, Basisobj.Rm);
                axang = [n theta];
            elseif isequal(size(rotm), [4, 1])
                axang = rotm; % Assuming it's already in axis-angle format
            else
                error('Invalid rotation matrix format.');
            end
            if isequal(t, [0, 0, 0])
                Umat = rotation(Basisobj.BasisL, axang, rightorleft);
            else
                % Handle translation here (if needed)
            end
        end
    end
end
