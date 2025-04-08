classdef Braid < matlab.mixin.CustomDisplay
    % Braid class for representing and processing braids
    %   This class encapsulates various properties and methods for working
    %   with braid words, their cycle decompositions, permutations, and
    %   generators. It also includes visualization and numerical methods.

    properties
        BraidWord;         % The braid word as a string
        Nstrings;          % Number of strings in the braid
        Crossings;         % Number of crossings in the braid
        NLink;             % Number of links in the braid
    end

    properties
        Generator;         % Generator representation of the braid
        BraidNum;          % Numerical representation of the braid word
        ChP;               % Placeholder for some braiding-related property (need more context)
        Kesi_k;            % Placeholder for another property
        Hsym;              % Symmetry information
        HTBkit;            % Placeholder for additional kit information
    end

    properties
        NaiveSeperateBands;  % Separated bands (likely related to eigenstates or physical modes)
        NaiveEIGENCAR;       % Eigenvalue/Car (perhaps eigenvectors or eigenvalue data)
        PlotSeqL;            % Plot sequence for braid representation
        PlotSeqLMirror;      % Mirrored plot sequence for braid representation
    end

    properties
        Permutation;            % Permutation associated with the braid
        CycleDecompositionStr;  % String representation of the cycle decomposition
        CycleDecomposition;     % The cycle decomposition structure
        CycleLift;              % Lifted cycle representation
    end

    properties % Data
        D_C = cell(1);      % Placeholder data cell
        kD_C = cell(1);     % Placeholder data cell
        ColorD_C = cell(1); % Placeholder for colored data
        Fnk = cell(1);      % Placeholder function/data
        Fjnk;               % Placeholder function/data
        kcross;             % Placeholder for crossing information
        Fkcross;            % Placeholder for crossing function
        kp;                 % Placeholder for momentum-related data
        Gkp;                % Placeholder for some Gkp data
        Gnk = cell(1);      % Placeholder for Gnk data
        Gjnk;               % Placeholder for Gjnk data
    end

    %% Protected Methods
    methods (Access = protected)
        % Define which properties should be displayed
        function propgrp = getPropertyGroups(~)
            proplist = {'BraidWord', 'Nstrings', 'Crossings', 'NLink', 'CycleDecompositionStr'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end

    %% Constructor Method
    methods
        function BraidObj = Braid(BraidWord, Nstrings)
            % Braid Constructor
            %   This constructor initializes a Braid object by processing
            %   the given braid word and the number of strings.

            % Convert braid word to numerical form
            DoubleBraidWord_col = Braid.BraidWord2Num(BraidWord);

            if nargin < 2
                % If Nstrings is not provided, determine it from the braid word
                Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
            end

            % Calculate the number of crossings from the braid word
            Crossings = length(DoubleBraidWord_col);

            % Generate braid-related data
            [Generator, GeneratorStr] = Braid.Nstrings2Generator(Nstrings);
            Permutation = Braid.PermutationBraidNum(DoubleBraidWord_col);
            [CycleDecomposition, CyclePresentationStr, CycleLift] = Braid.Cauchy2Cycle(Permutation);

            % Assign values to object properties
            BraidObj.BraidWord = BraidWord;
            BraidObj.BraidNum = DoubleBraidWord_col;
            BraidObj.Nstrings = Nstrings;
            BraidObj.Crossings = Crossings;
            BraidObj.Generator = Generator;
            BraidObj.Permutation = Permutation;
            BraidObj.CycleDecomposition = CycleDecomposition;
            BraidObj.CycleLift = CycleLift;
            BraidObj.NLink = length(CycleDecomposition);
            BraidObj.CycleDecompositionStr = CyclePresentationStr;

            % Generate initial bands and sequences
            [BraidObj.NaiveSeperateBands, BraidObj.NaiveEIGENCAR, BraidObj.PlotSeqL, BraidObj.PlotSeqLMirror] = NaiveBands(BraidObj);
        end
    end

    %% Methods
    methods
        function [Sign, k_index] = CheckCrossSign(BraidObj, kpoint)
            % CheckCrossSign: Determines the sign of the braid at a given k-point
            %   This method checks the sign of the braid at a specific k-point
            %   by calculating the index from the k-point and checking the
            %   corresponding sign of the braid number.

            BraidNum = BraidObj.BraidNum;
            Cross = length(BraidNum);
            BraidNumsign = sign(BraidNum);  % Sign of braid numbers
            k_normal = kpoint ./ (2 * pi) * Cross;  % Normalize k-point to braid scale
            k_index = ceil(k_normal);  % Find the index of the corresponding braid number
            Sign = BraidNumsign(k_index);  % Return the sign at that index
        end
    end

    %% Parent Function
    methods
        % Generate D_C, kD_C, and ColorD_C for the braid
        function BraidObj = D_C_gen(BraidObj)
            % Normalize the eigenvalues from NaiveEIGENCAR for trigonometric interpolation
            NormalPbands = BraidObj.NaiveEIGENCAR;
            Width = (max(NormalPbands, [], 'all') - min(NormalPbands, [], 'all')) / 2;
            Middle = (max(NormalPbands, [], 'all') + min(NormalPbands, [], 'all')) / 2;
            NormalPbands = -1 * (NormalPbands - Middle) / Width;

            % Initialize cycle decomposition and lift
            Lk = BraidObj.CycleDecomposition;
            CycleLift = BraidObj.CycleLift;
            Nsn = length(Lk);  % Number of strings

            % Initialize D_C as the normalized eigenvalues
            D_Cj = NormalPbands;
            D_C = cell(1, Nsn);
            ColorD_C = cell(1, Nsn);
            kD_C = cell(1, Nsn);

            % Loop over each cycle to build D_C and kD_C
            for n = 1:Nsn
                CycleTmpLift = CycleLift{n};
                s = length(CycleTmpLift);
                D_C{n} = D_Cj(CycleTmpLift(1), :);  % Start with the first element
                ColorD_C{n} = CycleTmpLift(1);  % Store the color

                for j = 2:s
                    D_C{n} = [D_C{n}, D_Cj(CycleTmpLift(j), 2:end)];  % Add subsequent elements
                    ColorD_C{n} = [ColorD_C{n}, CycleTmpLift(j)];  % Store color for each element
                end

                % Generate kD_C for Fourier interpolation
                kD_C{n} = linspace(0, 2*pi, size(D_C{n}, 2));
            end

            % Store results in BraidObj
            BraidObj.D_C = D_C;
            BraidObj.kD_C = kD_C;
            BraidObj.ColorD_C = ColorD_C;
        end

        % Generate the Fourier coefficients (Fnk) for the braid
        function [BraidObj, am_n] = FnkGen(BraidObj)
            % Check if the Fourier coefficients are already computed
            if isempty(BraidObj.Fjnk)
                BraidObj = D_C_gen(BraidObj);  % Generate D_C data if not already done
            end

            Lk = BraidObj.CycleDecomposition;
            Cross = BraidObj.Crossings;
            am_n = cell(1, numel(BraidObj.D_C));

            % Loop through each cycle for Fourier transformation
            for n = 1:numel(BraidObj.D_C)
                ln = length(Lk{n});
                Crossln = Cross * ln;
                kD_Cn = BraidObj.kD_C{n}(1:end-1);  % Remove last element
                D_Cn = BraidObj.D_C{n}(1:end-1);  % Remove last element

                % Adjust for odd or even cycle length
                if mod(Crossln, 2) == 1
                    mL = (-Crossln / 2 + 1 / 2) : (Crossln / 2 - 1 / 2);
                else
                    mL = (-Crossln / 2 + 1) : (Crossln / 2 - 1);
                    m = Crossln / 2;
                    expL = exp(-1i * kD_Cn * m);
                    a_Crossln_2 = 1 / Crossln * sum(D_Cn .* expL);
                end

                % Initialize and compute the Fourier coefficients
                am_n{n} = zeros(1, length(mL));
                for count = 1:length(mL)
                    expL = exp(-1i * kD_Cn * mL(count));
                    am_n{n}(count) = 1 / Crossln * sum(D_Cn .* expL);
                end

                % Create the Fourier transform expression
                syms k real;
                expkL = exp(1i * mL * k);
                if mod(Crossln, 2) == 1
                    Fnk{n} = simplify(rewrite(sum(expkL .* am_n{n}), 'sincos'));
                else
                    Fnk{n} = simplify(rewrite(sum([expkL .* am_n{n}, a_Crossln_2 * cos(Crossln / 2 * k)]), 'sincos'));
                end

                % Compute values at specific k-points
                kjn = (k + 2 * pi * (0:(ln - 1))) / ln;
                for ikj = 1:numel(kjn)
                    Fjnk{n}{ikj} = subs(Fnk{n}, k, kjn(ikj));
                end
            end

            % Store the Fourier transforms in BraidObj
            BraidObj.Fnk = Fnk;
            BraidObj.Fjnk = Fjnk;
        end

        % Generate k-point data for the braid's crossings and Fourier coefficients
        function BraidObj = kpGen(BraidObj, kn)
            % Optional parameter kn with default value 100
            if nargin < 2
                kn = 100;
            end

            % Generate Fnk data if not already computed
            if isempty(BraidObj.kcross)
                BraidObj = FnkGen(BraidObj);
            end

            % Initialize variables and k-point calculations
            Lk = BraidObj.CycleDecomposition;
            CycleLift = BraidObj.CycleLift;
            Falljn_label = [];
            count = 0;

            % Create a list of Fjnk data for the braid
            Falljn = {};
            for n = 1:numel(BraidObj.Fjnk)
                ln = length(BraidObj.Fjnk{n});
                for j = 1:ln
                    count = count + 1;
                    Falljn{count} = BraidObj.Fjnk{n}{j};
                    Falljn_label = [Falljn_label; [n, j, ln]];
                end
            end

            % Initialize arrays for k-crossings and the Fourier coefficients
            kcross = cell(1, numel(Falljn));
            Fkcross = cell(1, numel(Falljn));
            kp = cell(1, numel(Falljn));
            Gkp = cell(1, numel(Falljn));
            nF = count;

            % Loop over pairs of Fjnk and find crossing locations
            for i = 1:nF
                TargetFjnk_left = Falljn{i};
                n = Falljn_label(i, 1);
                CycleLiftn = CycleLift{n};
                jforstring = Falljn_label(i, 2);
                iString = CycleLiftn(jforstring);
                n_using = n;
                ln_using = Falljn_label(i, 3);
                j_using = jforstring;

                % Compare with every other Fjnk to find crossings
                for j = 1:nF
                    if i == j
                        continue;
                    end
                    TargetFjnk_right = Falljn{j};
                    n = Falljn_label(j, 1);
                    CycleLiftn = CycleLift{n};
                    jforstring = Falljn_label(j, 2);
                    jString = CycleLiftn(jforstring);

                    % Solve for crossing locations
                    EQtmpZero = simplify(TargetFjnk_left - TargetFjnk_right);
                    EQtmp = EQtmpZero == 0;
                    EQtmpF = matlabFunction(EQtmpZero);
                    EQtmpL = EQtmpF(linspace(0, 2*pi, kn));
                    EQtmpsignL = sign(EQtmpL);
                    EQtmpsigndiffL = diff(EQtmpsignL);
                    CrossLocation = find(EQtmpsigndiffL);

                    if ~isempty(CrossLocation)
                        for d = CrossLocation
                            syms k real;
                            EQtmp2 = (k > sym(linspace(0, 2*pi, kn)*(d-1)));
                            EQtmp3 = (k < sym(linspace(0, 2*pi, kn)*(d+1)));
                            kcrosstmp = solve([EQtmp, EQtmp2, EQtmp3], k);

                            % Handle possible solver failures
                            if isempty(kcrosstmp)
                                kcrosstmp = vpasolve(EQtmp, k, [linspace(0, 2*pi, kn)*(d-1), linspace(0, 2*pi, kn)*(d+1)]);
                            end

                            Fkcrosstmp = subs(TargetFjnk_left, k, kcrosstmp);
                            [TmpSign, k_index] = BraidObj.CheckCrossSign(kcrosstmp);
                            StringSeq = BraidObj.NaiveEIGENCAR([iString, jString]);

                            kcross{d} = kcrosstmp;
                            Fkcross{d} = Fkcrosstmp;
                            kp{d} = StringSeq;
                            Gkp{d} = TmpSign;
                        end
                    end
                end
            end

            % Store k-point data in BraidObj
            BraidObj.kcross = kcross;
            BraidObj.Fkcross = Fkcross;
            BraidObj.kp = kp;
            BraidObj.Gkp = Gkp;
        end

        % Generate kpset and Gkpset for the braid
        function BraidObj = GkData_gen(BraidObj)
            % kpset holds the k-crossing points, which are modified to generate kpsetforGk
            kpset = BraidObj.kcross;
            kpsetforGk = kpset;  % Initialize kpsetforGk
            syms k real;

            % Loop over each crossing point set
            for n = 1:numel(kpset)
                ln = length(BraidObj.CycleDecomposition{n});
                tmpkset = [];

                % Generate a set of k-values (adjusted by multiples of 2*pi)
                for j = 0:ln-1
                    tmpkset = [tmpkset, (kpset{n} + 2*pi*j) / ln];
                end

                % Check the sign at each crossing point
                for i = 1:numel(tmpkset)
                    Gkpset{n}(i) = BraidObj.CheckCrossSign(tmpkset(i));
                end

                kpsetforGk{n} = tmpkset;  % Update kpsetforGk with the new k-values
            end

            % Store the results in the Braid object
            BraidObj.kpsetforGk = kpsetforGk;
            BraidObj.Gkpset = Gkpset;
        end

        % Generate Gnk and Gjnk using the method from PHYSICAL REVIEW LETTERS 126, 010401 (2021)
        function [BraidObj, bm_n] = GnkGen2(BraidObj)
            % This method computes Gnk using trigonometric interpolation and an adjusted D_C
            % Check if Gkp is available, if not, generate it
            if isempty(BraidObj.Gkp)
                BraidObj = kpGen(BraidObj);
            end
            Lk = BraidObj.CycleDecomposition;

            % Loop through each cycle
            for n = 1:numel(BraidObj.Gkp)
                yk = BraidObj.Gkp{n};
                ln = length(Lk{n});
                N = length(yk);
                tkL = BraidObj.kp{n};  % k-point data for the cycle
                Odd = mod(N, 2) == 1;  % Check if the length is odd
                jL = 1:N;  % Index for yk values
                syms k real;
                expL = exp(1i * tkL);  % Exponential function for the k-points
                ChooseL = 1:N;

                % Handle odd or even N
                if Odd
                    K = (N - 1) / 2;
                    CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                    count = 0;

                    % Loop through each k-point for trigonometric interpolation
                    for tk_prime = tkL
                        count = count + 1;
                        ChooseLTmp = ChooseL;
                        ChooseLTmp(ChooseL == count) = [];  % Remove the current index
                        expLmodify = expL(ChooseLTmp);  % Adjust the exponential term
                        Cumprod(count) = fold(@times, (exp(1i * k) - expLmodify) ./ (exp(1i * tk_prime) - expLmodify));
                    end
                    Gnk{n} = simplify(sum(CoeffL .* Cumprod));  % Final expression for Gnk
                else
                    K = N / 2;
                    CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                    count = 0;

                    % Handle even N case
                    for tk_prime = tkL
                        count = count + 1;
                        ChooseLTmp = ChooseL;
                        ChooseLTmp(ChooseL == count) = [];
                        expLmodify = expL(ChooseLTmp);
                        expLmodify = [1, expLmodify];  % Adjust for the even case
                        Cumprod(count) = fold(@times, (exp(1i * k) - expLmodify) ./ (exp(1i * tk_prime) - expLmodify));
                    end
                    Gnk{n} = simplify(sum(CoeffL .* Cumprod));  % Final expression for Gnk
                end

                % Compute k-point values and substitute into the final Gnk expression
                kjn = (k + 2 * sym(pi) * (0:(ln - 1))) / ln;
                for ikj = 1:numel(kjn)
                    Gjnk{n}{ikj} = subs(Gnk{n}, k, kjn(ikj));  % Substitute k-values
                end
            end

            % Store the results in the Braid object
            BraidObj.Gnk = Gnk;
            BraidObj.Gjnk = Gjnk;
        end

        % Generate Gnk and Gjnk using the method from PHYSICAL REVIEW LETTERS 126, 010401 (2021)
        function [BraidObj, bm_n] = GnkGen(BraidObj)
            % This method generates Gnk and Gjnk using a trigonometric approach
            if isempty(BraidObj.Gkp)
                BraidObj = kpGen(BraidObj);  % Ensure Gkp is available
            end
            Lk = BraidObj.CycleDecomposition;

            % Loop through each cycle for Fourier analysis
            for n = 1:numel(BraidObj.Gkp)
                D_Cn = BraidObj.Gkp{n};
                ln = length(Lk{n});
                N = length(D_Cn);  % Number of elements in the cycle
                kD_Cn = BraidObj.kp{n};
                Odd = mod(N, 2) == 1;
                syms k real;

                % Handle odd and even cases for the Fourier series
                Crossln = N;
                Crossln_2 = N / 2;
                if Odd
                    mL = (-Crossln / 2 + 1 / 2):(Crossln / 2 - 1 / 2);
                    Width = length(mL);
                    expSymL = exp(1i * k .* mL);
                    expLFunction = matlabFunction(expSymL, 'Vars', k);  % Convert to function for efficiency
                else
                    mL = (-Crossln / 2 + 1):(Crossln / 2 - 1);
                    Width = length(mL) + 1;
                    expSymL = [exp(1i * k .* mL), cos(Crossln_2 * k)];
                    expLFunction = matlabFunction(expSymL, 'Vars', k);
                end

                % Create the matrix for the system of equations
                TheMat = zeros([N, Width], 'sym');
                for i = 1:N
                    TheMat(i, :) = expLFunction(kD_Cn(i));
                end
                TheMat = simplify(TheMat);  % Simplify the matrix

                % Check if the determinant is zero, which indicates a degenerate case
                if double(det(TheMat)) == 0
                    yk = BraidObj.Gkp{n};
                    ln = length(Lk{n});
                    tkL = BraidObj.kp{n};
                    jL = 1:N;
                    expL = exp(1i .* tkL);
                    ChooseL = 1:N;

                    % Handle degenerate case using trigonometric interpolation
                    if Odd
                        K = (N - 1) / 2;
                        CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                        count = 0;

                        for tk_prime = tkL
                            count = count + 1;
                            ChooseLTmp = ChooseL;
                            ChooseLTmp(ChooseL == count) = [];
                            expLmodify = expL(ChooseLTmp);
                            Cumprod(count) = fold(@times, (exp(1i * k) - expLmodify) ./ (exp(1i * tk_prime) - expLmodify));
                        end
                        Gnk{n} = simplify(sum(CoeffL .* Cumprod));
                    else
                        K = N / 2;
                        CoeffL = yk(jL) .* exp((-1i * K * k) + (1i * K * tkL));
                        count = 0;

                        for tk_prime = tkL
                            count = count + 1;
                            ChooseLTmp = ChooseL;
                            ChooseLTmp(ChooseL == count) = [];
                            expLmodify = expL(ChooseLTmp);
                            expLmodify = [1, expLmodify];
                            Cumprod(count) = fold(@times, (exp(1i * k) - expLmodify) ./ (exp(1i * tk_prime) - expLmodify));
                        end
                        Gnk{n} = simplify(sum(CoeffL .* Cumprod));
                    end
                else
                    % Solve the system using Cramer's rule if the matrix is non-singular
                    TheSolve = (TheMat) \ D_Cn.';
                    Gnk{n} = simplify(sum((TheSolve.').* expSymL));
                end

                % Substitute k values into the expression for Gjnk
                kjn = (k + 2 * sym(pi) * (0:(ln - 1))) / ln;
                for ikj = 1:numel(kjn)
                    Gjnk{n}{ikj} = subs(Gnk{n}, k, kjn(ikj));
                end
            end

            % Store the results in the Braid object
            BraidObj.Gnk = Gnk;
            BraidObj.Gjnk = Gjnk;
        end

        % Generate characteristic polynomial (ChP) for the braid
        function BraidObj = ChPGen(BraidObj)
            syms lambda k real;
            Flambdak = sym(1);  % Initialize characteristic polynomial

            % Generate Gnk if not already generated
            if isempty(BraidObj.Gjnk)
                BraidObj = BraidObj.GnkGen();
            end

            % Fjnk and Gjnk contain the computed values for the braid
            FjnCell = BraidObj.Fjnk;
            GjnCell = BraidObj.Gjnk;

            % Loop through each link and each cycle to build the characteristic polynomial
            for n = 1:BraidObj.NLink
                ln = length(BraidObj.CycleDecomposition{n});
                for j = 1:ln
                    Flambdak = Flambdak * (lambda - FjnCell{n}{j} - 1i * GjnCell{n}{j});
                end
            end

            % Simplify the characteristic polynomial
            Flambdak = simplify(Flambdak);
            BraidObj.ChP = Flambdak;

            % Compute the Laurent series and store coefficients
            LaurentSeries = series(BraidObj.ChP, lambda, 'Order', BraidObj.Nstrings);
            BraidObj.Kesi_k = coeffs(LaurentSeries, lambda);
        end

        % Generate Hsym matrix (symmetry matrix) from characteristic polynomial
        function BraidObj = HsymGen(BraidObj)
            syms k k_x real;

            % Generate ChP if not already generated
            if isempty(BraidObj.ChP)
                BraidObj = BraidObj.ChPGen();
            end

            % Construct the Hsym matrix (a Toeplitz-like structure)
            Hsym = sym(diag(ones(1, BraidObj.Nstrings - 1), -1));
            Hsym(1, :) = -flip(BraidObj.Kesi_k);  % Adjust the first row based on Kesi_k

            % Substitute k_x for k in the Hsym matrix (for k-space representation)
            BraidObj.Hsym = subs(Hsym, k, k_x);
        end

        % Generate Tight-Binding Kit (TBkit) for the braid
        function BraidObj = TBkitGen(BraidObj)
            syms k k_x real;

            % Initialize the Tight-Binding Kit object
            HTBkit = TBkit();
            HTBkit.Basis_num = BraidObj.Nstrings;
            HTBkit.Rm = [1];
            HTBkit.Dim = 1;
            HTBkit.Hsym = BraidObj.Hsym;  % Set the symmetry matrix from Hsym

            % Define k-points for the tight-binding model (in Cartesian coordinates)
            kpoints_cart = [
                0;
                pi;
                pi;
                2*pi;
                ];

            % Convert to fractional k-points (dividing by 2*pi)
            kpoints_frac = kpoints_cart / (2*pi);

            % Generate the k-path for the tight-binding kit with 61 k-points
            [HTBkit.klist_cart, HTBkit.klist_frac, HTBkit.klist_l, HTBkit.kpoints_l, ~] = ...
                TBkit.kpathgen(kpoints_frac, 61, HTBkit.Gk, 'Dim', 1);

            % Assign k-point names for the labels
            HTBkit.kpoints_name = ["0", "\pi", "2\pi"];

            % Store the generated HTBkit object in the Braid object
            BraidObj.HTBkit = HTBkit;
        end
    end
    %% DataSet
    methods
        function [NaiveSeperateBands, NaiveEIGENCAR, PlotSeqL, PlotSeqLVertical] = NaiveBands(BraidObj)
            % Initialize variables
            NaiveSeperateBands = cell(1, BraidObj.Crossings); % Cell array to store band separations
            PermutationList = abs(BraidObj.BraidNum); % Absolute value of braid numbers
            BraidList = BraidObj.BraidNum; % List of braid numbers
            Nbands = BraidObj.Nstrings; % Number of bands (strings)
            Ncross = BraidObj.Crossings; % Number of crossings

            % Default indices for band separation
            DefaultL = (1:Nbands).';

            % Initialize plot sequence matrices
            PlotSeqL = repmat(DefaultL, [1, Ncross]);
            PlotSeqLVertical = PlotSeqL; % For vertical plotting
            NaiveEIGENCAR = zeros(Nbands, Ncross + 1); % Store eigenvalue data
            NaiveEIGENCAR(:, 1) = DefaultL; % Initialize the first column with DefaultL

            % Loop through each crossing
            for i = 1:Ncross
                % Set up the first column of NaiveSeperateBands for this crossing
                if i > 1
                    NaiveSeperateBands{i}(:, 1) = NaiveSeperateBands{i - 1}(:, 2);
                else
                    NaiveSeperateBands{i}(:, 1) = DefaultL;
                end

                % Initialize flags for positive and negative band indices
                PlusN = -1;
                MinusN = -1;

                % Loop through each band and update band positions based on permutations
                for n = 1:Nbands
                    % Find bands that need to be incremented or decremented based on the permutation list
                    if NaiveSeperateBands{i}(n, 1) == PermutationList(i)
                        NaiveSeperateBands{i}(n, 2) = NaiveSeperateBands{i}(n, 1) + 1;
                        PlusN = n; % Remember the index for the positive band
                    elseif NaiveSeperateBands{i}(n, 1) == PermutationList(i) + 1
                        NaiveSeperateBands{i}(n, 2) = NaiveSeperateBands{i}(n, 1) - 1;
                        MinusN = n; % Remember the index for the negative band
                    else
                        NaiveSeperateBands{i}(n, 2) = NaiveSeperateBands{i}(n, 1);
                    end
                end

                % Swap the positions of bands if necessary based on braid direction
                if PlusN > 0 && MinusN > 0
                    if BraidList(i) > 0 && PlusN > MinusN
                        % Do nothing if braid is positive and the order is correct
                    elseif BraidList(i) < 0 && PlusN < MinusN
                        % Do nothing if braid is negative and the order is correct
                    else
                        % Swap the positions of PlusN and MinusN
                        PlotSeqL([MinusN, PlusN], i) = PlotSeqL([PlusN, MinusN], i);
                        PlotSeqLVertical([MinusN, PlusN], i) = PlotSeqLVertical([PlusN, MinusN], i);
                    end
                end

                % Store the updated bands into the NaiveEIGENCAR array
                NaiveEIGENCAR(:, i + 1) = NaiveSeperateBands{i}(:, 2);
            end
        end
    end
    %% plot
    methods
        function ax = Fjnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = [1/2,1/2,1/2];
                opts.vertical = false;
                opts.kn = 100;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            if isempty(BraidObj.Fkcross)
                BraidObj = BraidObj.kpGen();
            end
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            Fjnk = BraidObj.Fjnk;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*50;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    for j = 1:length(ColorD_C{n})
                        FjnkFunction = matlabFunction(Fjnk{n}{j},'Vars',sym('k'));
                        Ei = ColorD_C{n}(j);
                        KPs = [0,2*pi];
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FjnkFunction(KL);
                        plot(ax,double(FL),double(KL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['F_',num2str(n),'^',num2str(j),'(k)']);
                        SelectL = SelectL + Ndivide -1;
                    end
                    scatter(ax,...
                        double(BraidObj.Fkcross{n}),...
                        double(BraidObj.kcross{n}),...
                        options.MarkerSize,...
                        options.MarkerFaceColor,...
                        'filled','DisplayName','Cross');
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    for j = 1:length(ColorD_C{n})
                        FjnkFunction = matlabFunction(Fjnk{n}{j},'Vars',sym('k'));
                        Ei = ColorD_C{n}(j);
                        KPs = [0,2*pi];
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FjnkFunction(KL);
                        plot(ax,double(KL),double(FL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['F_',num2str(n),'^',num2str(j),'(k)']);
                        %fplot(ax,[0,2*pi],options.LineSpec,"LineWidth",options.LineWidth);
                        SelectL = SelectL + Ndivide -1;
                    end
                    scatter(ax,double(BraidObj.kcross{n}),...
                        double(BraidObj.Fkcross{n}),...
                        options.MarkerSize,...
                        options.MarkerFaceColor,...
                        'filled','DisplayName','Cross');
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = Fnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
                opts.kn = 100;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            Fnk = BraidObj.Fnk;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    FnkFunction = matlabFunction(Fnk{n},'Vars',sym('k'));
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        KPs = kD_C{n}(SelectL);
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FnkFunction(KL);
                        plot(ax,double(FL),double(KL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['EFit_',num2str(Ei)]);
                        plot(ax,double(D_C{n}(SelectL)),double(KPs),'o',...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    FnkFunction = matlabFunction(Fnk{n},'Vars',sym('k'));
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        KPs = kD_C{n}(SelectL);
                        KL = linspace(KPs(1),KPs(end),opts.kn);
                        FL = FnkFunction(KL);
                        plot(ax,double(KL),double(FL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'DisplayName',['EFit_',num2str(Ei)]);
                        plot(ax,double(KPs),D_C{n}(SelectL),'o',...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        %fplot(ax,[0,2*pi],options.LineSpec,"LineWidth",options.LineWidth);
                        SelectL = SelectL + Ndivide -1;
                    end
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = Gnk_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 5;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'r';
                opts.vertical = false;
                opts.kn = 361;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            Gkp = BraidObj.Gkp;
            kp = BraidObj.kp;
            Gnk = BraidObj.Gnk;
            [~,Ax] = Figs(length(BraidObj.Gkp),1);
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    ticks = double(kp{n});
                    ticklabel = repmat("",[1 length(kp{n})]);
                    for i = 1:length(kp{n})
                        ticklabel(i) = latex(kp{n}(i));
                    end
                    GnkFunction = matlabFunction(Gnk{n},'Vars',sym('k'));
                    KPs = [0,2*pi];
                    KL = linspace(KPs(1),KPs(end),opts.kn);
                    GL = GnkFunction(KL);
                    plot(ax,double(GL),double(KL),options.LineSpec,...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerEdgeColor,...
                        'DisplayName',['G_',num2str(n),'(k)']);
                    plot(ax,double(Gkp{n}),double(kp{n}),'o',...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerFaceColor ,...
                        'MarkerSize',options.MarkerSize,...
                        'MarkerEdgeColor',options.MarkerEdgeColor,...
                        'MarkerFaceColor',options.MarkerFaceColor ,...
                        'DisplayName',['G_',num2str(n),'(k_p)']);
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    %xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    ticks = double(kp{n});
                    ticklabel = repmat("",[1 length(kp{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kp{n})
                        ticklabel(i) = "$"+latex(kp{n}(i))+"$";
                    end
                    GnkFunction = matlabFunction(Gnk{n},'Vars',sym('k'));
                    KPs = [0,2*pi];
                    KL = linspace(KPs(1),KPs(end),opts.kn);
                    GL = GnkFunction(KL);
                    plot(ax,double(KL),double(GL),options.LineSpec,...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerEdgeColor,...
                        'DisplayName',['G_',num2str(n),'(k)']);
                    plot(ax,double(kp{n}),double(Gkp{n}),'o',...
                        'LineWidth',options.LineWidth,...
                        'Color',options.MarkerFaceColor ,...
                        'MarkerSize',options.MarkerSize,...
                        'MarkerEdgeColor',options.MarkerEdgeColor,...
                        'MarkerFaceColor',options.MarkerFaceColor ,...
                        'DisplayName',['G_',num2str(n),'(k_p)']);
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    %ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = D_C_Plot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                %options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-o';
                options.LineWidth = 10;
                options.MarkerSize = -1;
                options.MarkerEdgeColor = [1/2,1/2,1/2];
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
            end
            Nbands=BraidObj.Nstrings;
            %BraidObj = D_C_gen(BraidObj);
            D_C = BraidObj.D_C;
            kD_C = BraidObj.kD_C;
            ColorD_C = BraidObj.ColorD_C;
            [~,Ax] = Figs(length(BraidObj.D_C),1);
            CycleLift = BraidObj.CycleLift;
            %ax = options.ax;
            %hold(ax,'on');
            if options.MarkerSize <0
                options.MarkerSize  = options.LineWidth*3;
            end
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for n = 1:Nbands
                    colorL(n,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);

            if opts.vertical
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = size(D_C{n},2)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    %CycleLiftn = CycleLift{n};
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = latex(kD_C{n}(i));
                    end
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        plot(ax,D_C{n}(SelectL),kD_C{n}(SelectL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    ylabel(ax,'k');
                    yticks(ax,ticks);
                    set(ax,'TickLabelInterpreter','latex');
                    set(ax,'FontName','Times');
                    yticklabels(ax,ticklabel);
                    ylim(ax,[0-pi/10,2*pi+pi/10]);
                    xlim(ax,[-1.2,1.2]);
                    xlabel(ax,'x');
                    xticks(ax,[-1,0,1]);
                end
            else
                for n = 1:length(Ax)
                    ax = Ax(n);
                    Ndivide = 1+ (size(D_C{n},2)-1)/length(ColorD_C{n});
                    SelectL = 1:Ndivide;
                    ticks = double(kD_C{n});
                    ticklabel = repmat("",[1 length(kD_C{n})]);
                    set(ax,'FontName','Times');
                    for i = 1:length(kD_C{n})
                        ticklabel(i) = "$"+latex(kD_C{n}(i))+"$";
                    end
                    for j = 1:length(ColorD_C{n})
                        Ei = ColorD_C{n}(j);
                        plot(ax,double(kD_C{n}(SelectL)),D_C{n}(SelectL),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',ModifycolorL(Ei,:),...
                            'DisplayName',['E_',num2str(Ei)]);
                        SelectL = SelectL + Ndivide -1;
                    end
                    set(ax,'TickLabelInterpreter','latex');
                    ylabel(ax,'y');
                    yticks(ax,[-1,0,1]);
                    xlabel(ax,'k');
                    xticks(ax,ticks);
                    xlim(ax,[0-pi/10,2*pi+pi/10]);
                    ylim(ax,[-1.2,1.2]);
                    xticklabels(ax,ticklabel);
                end
            end

            %axis(ax,'off');
        end
        function ax = LinePlot(BraidObj,options,opts)
            arguments
                BraidObj Braid
                options.ax = gca;
                options.Color = @parula;
                options.LineSpec = '-';
                options.LineWidth = 10;
                options.MarkerSize = 3;
                options.MarkerEdgeColor = 'none';
                options.MarkerFaceColor = 'none';
                opts.vertical = false;
            end
            Nbands=BraidObj.Nstrings;
            ax = options.ax;
            hold(ax,'on');
            if isstring(options.Color)
                if isrow
                    colorL = options.Color.';
                end
            elseif isnumeric(options.Color)
                x = linspace(0,1,size(options.Color,1));
                xq = linspace(0,1,Nbands);
                colorL = [...
                    interp1(x,options.Color(:,1),xq).',...
                    interp1(x,options.Color(:,2),xq).',...
                    interp1(x,options.Color(:,3),xq).',...
                    ];
            elseif isa(options.Color,'function_handle')
                colorL = options.Color(Nbands);
            else
                for i = 1:Nbands
                    colorL(i,:) = [rand,rand,rand];
                end
            end
            HSV = rgb2hsv(colorL);
            HSV(:,2) = HSV(:,2);
            HSV(:,3) = HSV(:,3)-0.2;
            ModifycolorL = hsv2rgb(HSV);
            PlotSeqList = BraidObj.PlotSeqL;
            EIGENCARCELL = BraidObj.NaiveSeperateBands;
            EIGENCAR = BraidObj.NaiveEIGENCAR;
            Ncross = BraidObj.Crossings;
            if opts.vertical
                PlotSeqList = BraidObj.PlotSeqL;
                for i = 1:Ncross
                    for n = 1:BraidObj.Nstrings

                        Ei = PlotSeqList(n,i);
                        plot(ax,EIGENCARCELL{i}(Ei,:),[i-1 i],options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,EIGENCAR(n,:),0:Ncross,...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            else
                for i = 1:Ncross
                    EIGENCAR(:,i+1) = EIGENCARCELL{i}(:,2);
                    for n = 1:BraidObj.Nstrings

                        Ei = PlotSeqList(n,i);
                        plot(ax,[i-1 i],-EIGENCARCELL{i}(Ei,:),options.LineSpec,...
                            'LineWidth',options.LineWidth,...
                            'Color',colorL(Ei,:),...
                            'MarkerSize',options.MarkerSize,...
                            'MarkerEdgeColor',options.MarkerEdgeColor,...
                            'MarkerFaceColor',options.MarkerFaceColor,...
                            'DisplayName',['E_',num2str(Ei)]);
                    end
                end
                for n = 1:Nbands
                    scatter(ax,0:Ncross,-EIGENCAR(n,:),...
                        ones(1,Ncross+1)*options.LineWidth*20,ModifycolorL(n,:),"filled",'DisplayName',['String_',num2str(n)]);
                end
            end
            ylabel(ax,'');
            yticks(ax,[]);
            xlabel(ax,'');
            xticks(ax,[]);
            axis(ax,'off');
        end
    end
    %% Construct Tools
    methods(Static)
        function BraidNum = BraidWord2Num(BraidWord)
            CharBraidWord = char(BraidWord);
            CharBraidWord_col = (CharBraidWord.');
            BraidNum = double(CharBraidWord_col);
            for i = 1:numel(BraidNum)
                if BraidNum(i) < 97
                    BraidNum(i) =  BraidNum(i) - 64;
                else
                    BraidNum(i) = -BraidNum(i) + 96;
                end

            end

        end
        function Nstrings = NumBraidWord2Nstrings(DoubleBraidWord_col)
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
        end
        function [Generator,GeneratorStr] = Nstrings2Generator(Nstrings)
            TmpList =  string((1:Nstrings-1).');
            GeneratorStr = "sigma_" + TmpList;
            Generator = str2sym(GeneratorStr);
        end
        function Permutation = PermutationBraidNum(DoubleBraidWord_col)
            PermutationList = abs(DoubleBraidWord_col);
            Nstrings = max(abs(DoubleBraidWord_col)) + 1;
            Permutation(1,:) = 1:Nstrings;
            for iString = Permutation(1,:)
                Permutation(2,iString) = Braid.FinalPosition(iString,PermutationList);
            end
        end
        function oString = FinalPosition(iString,PermutationList)
            tmpString =  iString;
            for i = 1:numel(PermutationList)
                if tmpString == PermutationList(i)
                    tmpString = tmpString + 1;
                elseif tmpString == PermutationList(i) +1
                    tmpString = tmpString - 1;
                else

                end
            end
            oString = tmpString;
        end
        function [CycleDecomposition,CycleDecompositionStr,CycleLift] = Cauchy2Cycle(Permutation)
            %CyclePresentation{1} = 1;
            PermutationStore = Permutation(1,:);
            PermutationFunction = Permutation(2,:);
            Ndivide = 0;
            while ~isempty(PermutationStore)
                TriceElement = PermutationStore(1);
                Converge = false;
                RefTriceElement = TriceElement;
                Ndivide = Ndivide + 1;
                count = 0;
                while ~Converge
                    count = count + 1;
                    CycleDecomposition{Ndivide}(count) = TriceElement;
                    PermutationStore(PermutationStore == TriceElement) = [];
                    TriceElement = PermutationFunction(TriceElement);
                    Converge = TriceElement == RefTriceElement;
                end
            end
            CycleDecompositionStr = '';
            for i = 1:numel(CycleDecomposition)
                CycleDecompositionStr = [CycleDecompositionStr,'('];
                CycleDecompositionStr = [CycleDecompositionStr,num2str(CycleDecomposition{i})];
                CycleDecompositionStr = [CycleDecompositionStr,')'];
                CycleLift{i} = CycleDecomposition{i};
            end

        end
        function Flambdak = f_lambdak(FnCell,GnCell,knCell,Ncycle,ln)
            syms lambda;
            syms k real;
            Flambdak = sym(1);
            for n = 1:Ncycle
                for j = 0:ln-1
                    Flambdak = Flambdak*(lambda-FnCell{n}(knCell{j+1})-1i*GnCell{n}(knCell{j+1}));
                end
            end

        end
    end
    methods(Static) % script
        function Permutation = PermutationBraidWord(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Permutation = Braid.PermutationBraidNum(DoubleBraidWord_col);
        end
        function Nstrings = BraidWord2Nstrings(BraidWord)
            DoubleBraidWord_col = Braid.BraidWord2Nstrings(BraidWord);
            Nstrings = Braid.NumBraidWord2Nstrings(DoubleBraidWord_col);
        end
    end
end