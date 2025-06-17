classdef Spin < HollowKnight
    properties(Access = private, Dependent)
        s
        ms
    end
    properties(Hidden)
        J
        Jz
    end
    properties
        parity = 1;
        orientation = [0 0 1];
    end
    properties(Dependent, Hidden)
        hollow;
        s2_bar;
        EigenVectors;
        EigenVector;
    end
    methods
        function SpinObj = Spin(S,Sz,coe,options)
            arguments
                S {mustBeHalfInteger(S)} = 0;
                Sz = [];
                coe = 1;
                options.orientation = [0 0 1];
                options.parity = 1;
                %options.coe = 1;
            end
            optionsCell = namedargs2cell(options);
            % call father : HollowKnight
            SpinObj = SpinObj@HollowKnight();
            %
            if isempty(S)
                SpinObj = Spin.empty([0 0]);
                return;
            end
            if isempty(Sz)
                %SpinObj.J = S;
                %SpinObj.orientation = options.orientation;
                Sz = (S:-1:-S).';
            end
            % use S = -1;Sz = nan;coe =nan; any of them repersent a empty
            % function

            % single
            SpinObj(1,1).J = S(1,1);
            SpinObj(1,1).orientation = options.orientation;
            if isinteger(SpinObj(1,1).J)
                SpinObj(1,1).parity = (-1)^SpinObj(1,1).J;
            else
                SpinObj(1,1).parity = options.parity;
            end
            SpinObj(1,1).Jz = Sz(1,1);
            SpinObj(1,1).coe = coe(1,1);

            % multi
            [nS_row,nS_col] = size(S);
            [nSz_row,nSz_col] = size(Sz);
            [ncoe_row,ncoe_col] = size(coe);
            nSpin = max([nS_row ,nSz_row, ncoe_row] );
            nSpin_col = max([nS_col ,nSz_col, ncoe_col] );
            if nSpin == 1 && nSpin_col == 1
                return;
            else
                SpinObj = repmat(SpinObj(1,1),[nSpin,nSpin_col]);
                [S,Sz,coe] = Spin.StandardInput(S,Sz,coe,'nrow',nSpin,'ncol',nSpin_col);
                for i = 1:numel(SpinObj)
                    SpinObj(i) = Spin(S(i),Sz(i),coe(i),optionsCell{:});
                end
            end
        end
    end
    %% overload
    methods
        A = uminus(A)
        C = minus(A, B, options)
        C = plus(A, B, options)
        C = innertimes(A, B)
        C = eq(A, B, options)
    end
    %% contract
    methods
        B = contractrow(A, options)
    end
    %% set
    methods
        SpinObj = setparity(SpinObj, paritymat)
    end
    %% get
    methods
        function hollow = get.hollow(SpinObj)
            if isnan(SpinObj.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
        function s = get.s(SpinObj)
            s = SpinObj.J;
        end
        function ms = get.ms(SpinObj)
            ms = SpinObj.Jz;
        end
        function s2_bar = get.s2_bar(SpinObj)
            s2_bar = SpinObj.J*(SpinObj.J+1);
        end
        %         function parity = get.parity(SpinObj)
        %             parity = (-1)^(SpinObj.J);
        %         end
        function EigenVectors = get.EigenVectors(SpinObj)
            %for i = 1:length(SpinObj)
            i=1;
            SzM = Sz(SpinObj(i),'full',true);
            SxM = Sx(SpinObj(i),'full',true);
            SyM = Sy(SpinObj(i),'full',true);
            [W,U] = eig(SxM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,1} = V/norm(V);
            [W,U] = eig(SyM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,2} = V/norm(V);
            [W,U] = eig(SzM);
            V= W(:,(diag(U) == SpinObj(i).Jz));
            EigenVectors{i,3} = V/norm(V);
            %end
        end
    end
    %% rotation
    methods
        U = rotation(SpinObj, rotm, rightorleft, options)
        U = rotation2(SpinObj, rotm, rightorleft, options)
    end
    
    methods
        Ak = rotateinner(A, abc, RightorLeft, immproper, conjugate, antisymmetry)
        Trmat = Tr(SpinObj)
        Invmat = ParityMat(SpinObj)
        WigerDmat = WignerD(SpinObj, abc, rightorleft)
        Matelement = WignerD_single(SpinObj1, SpinObj2, abc, rightorleft)
    end
    %% operator
    methods
        SpinObj = CG(Spinobj1, Spinobj2, options)
        SpinObj = TimeRerversal(SpinObj)
        SpinObj = Inversion(SpinObj)
        SpinObj = SzOper(Spinobj)
        SpinObj = SxOper(Spinobj)
        SpinObj = SyOper(Spinobj)
        SpinObj = SplusOper(Spinobj)
        SpinObj = SminusOper(Spinobj)
    end
    %% operator mat
    methods
        SzM = Sz(SpinObj,options)
        SyM = Sy(SpinObj,options)
        SxM = Sx(SpinObj,options)
        LzM = Lz(SpinObj,options)
        LxM = Lx(SpinObj,options)
        LyM = Ly(SpinObj,options)
        JzM = JZ(SpinObj,options)
        JxM = JX(SpinObj,options)
        JyM = JY(SpinObj,options)
        Splus = Splus(SpinObj,options)
        Sminus = Sminus(SpinObj,options)
        Lplus = Lplus(SpinObj,options)
        Lminus = Lminus(SpinObj,options)
        Jplus = Jplus(SpinObj,options)
        Jminus = Jminus(SpinObj,options)
    end
    %% Modify
    methods
        SpinObj = BasisGen(Spinobj,force)
        SpinObj = PerM(Spinobj,permutation)
    end
    %% pretty
    methods
        disp(SpinObj)
        Output = string(SpinObj)
        Output = pretty(SpinObj,options)
    end
    methods(Static)
    end
    methods(Static)
    end
end
