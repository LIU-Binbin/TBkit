classdef Y_l__m < HollowKnight
    properties
        l   double;
        m   double;
    end
    properties(Dependent)
        expression ;
        parity ;
    end
    properties(Dependent,Hidden)
        hollow;
    end
    properties(Hidden)
        n   = 0         ;
        symbolic = true ;
        explicit = false;
    end
    methods
        function YlmObj = Y_l__m(l,m,coe,n,TH,PHI,options)
            arguments
                l   = 1;%... Azimuthal quantum number
                m   = [];%... Magnetic quantum number
                coe = 1;
                n   = 0;%... Radius(n?)
                TH  = sym('theta','real');%... azimutal angles
                PHI = sym('phi','real');%... polar angles
                options.common = true;
                options.symbolic = false;
                options.normalization {mustBeMember(options.normalization,{'unnorm','norm','sch'})}= 'unnorm';
            end

            % call father : HollowKnight
            YlmObj = YlmObj@HollowKnight();

            % single
            YlmObj.symbolic = options.symbolic ;
            if ~options.common
                if isa(l,'sym') && isa(m,'sym')
                    switch options.normalization
                        case 'unnorm'
                            YlmObj.expression = sym('');
                        case 'norm'
                            YlmObj.expression = sym('');
                        case 'sch'
                    end
                else

                end
                if isa(TH,'sym')

                end
                if isa(PHI,'sym')
                end
            end
            if isempty(m)
                m = (l:-1:-l).';
            end
            YlmObj(1,1).l = l(1,1);
            YlmObj(1,1).m = m(1,1);
            YlmObj(1,1).coe = coe(1,1);
            YlmObj(1,1).n = n(1,1);

            % multi
            [nS_row,nS_col] = size(l);
            [nSz_row,nSz_col] = size(m);
            [ncoe_row,ncoe_col] = size(coe);
            [nn,nn_col] = size(n);
            nYlmObj = max([nS_row ,nSz_row, ncoe_row,nn]);
            nYlmObj_col = max([nS_col ,nSz_col, ncoe_col nn_col] );
            if nYlmObj == 1 && nYlmObj_col == 1
                return;
            else
                YlmObj = repmat(YlmObj(1,1),[nYlmObj nYlmObj_col]);
                [l,m,coe,n] = Y_l__m.StandardInput(l,m,coe,n,'nrow',nYlmObj,'ncol',nYlmObj_col);
                for i =1:numel(YlmObj)
                    YlmObj(i).l = l(i);
                    YlmObj(i).m = m(i);
                    YlmObj(i).n = n(i);
                    YlmObj(i).coe = coe(i);
                end
            end
        end
    end

    methods %get
        function hollow = get.hollow(YlmObj)
            if isnan(YlmObj.coe)
                hollow = true;
            else
                hollow = false;
            end
        end
        function parity = get.parity(YlmObj)
            if length(YlmObj) == 1
                parity = (-1)^YlmObj.l;
            else
                parity = YlmObj(1).parity;
                for i =2:length(YlmObj)
                    if parity ~=YlmObj(i).parity
                        parity = 0 ;
                        return;
                    end
                end
            end
        end
        function expression = get.expression(YlmObj)
            if YlmObj.l>4
                syms theta phi real;
                if YlmObj.symbolic
                    Symble = vasplib.SymbolicVarible('N',YlmObj.m,YlmObj.l);
                    expression = Symble...
                        *((-1)^YlmObj.m*(sqrt(1/2/sym(pi)))*exp(1i*YlmObj.m*phi));
                    expression = simplify(YlmObj.coe*expression);
                else
                    Symble = vasplib.SymbolicVarible('N',YlmObj.m,YlmObj.l);
                    expression = Symble...
                        *vpa((-1)^YlmObj.m*sqrt(1/2/pi)*exp(1i*YlmObj.m*phi));
                end
            else
                expression = explicitformula(YlmObj,'vpa',false);
                expression = simplify(YlmObj.coe*expression);
            end
        end
    end

    methods
        C = plus(A,B)
        A = umius(A)
        C = minus(A,B)
        C = innertimes(A,B)
        C = eq(A,B,options)
        C = mrdivide(A,B)
        C = lt(A,B,options)
        sum()
    end
    methods
        B = contractrow(A,options)
        B = LplusOper(A)
        B = LminusOper(A)
        B = LzOper(A)
    end
    methods
        Mat = Lplus(A)
        Mat = Lz(A)
        Mat = Lminus(A)
    end
    methods
        Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
    end
    methods
        disp(YlmObj,options)
        expression = pretty(YlmObj,options)
        Tesseralexpansion = Tesseral(Y_lmObj,options)
        str = string(YlmObj,options)
        expression = formula(YlmObj,options)
        expression = explicitformula(YlmObj,options)
    end
    methods(Static)
        SymExpr = d(j,m1,m2,seed)
        SymExpr = d_mm__j(j,m1,m2,seed)
        PlmExpr = Plm(l,m,seed,options)
        NlmExpr = Nlm(l,m,seed,options)
        YlmExpr = InnerY_l__m(l,m,seed1,seed2,options)
    end
    methods(Static)
        DFact = DFactorial(n)
        c = binomial(n, k)
        cg = CG(j1,m1,j2,m2,j,m)
        W = w3j(j1,j2,j3,m1,m2,m3)
        W = w6j(a,b,c,d,e,f)
        W = w9j(a,b,c,d,e,f,g,h,j)
        cg = ClebschGordan(j1,j2,j,m1,m2,m)
        W = Wigner3j( j123, m123 )
        W = Wigner3j_sym(j123, m123)
    end
    methods(Static)
        OutExpr = nlm2atomic(l,m,n,options)
        Outstr = l__m2str(l,m)
        Outstr = lm2str(l,m)
        [dupNames, dupNdxs] = getDuplicates(aList)
        Ind = IndexOfMultiples(A)
        T = isMultiple(A)
    end
    methods(Static)
        tri = del(a,b,c)
        tf = triangular_cond(a,b,c)
        sym_fact = sym_fact(n)
    end
end
