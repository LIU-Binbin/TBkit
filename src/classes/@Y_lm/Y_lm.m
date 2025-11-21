classdef Y_lm < Y_l__m
    %Y_lm  Real Spherical harmonic function.
    % rewrite this class; is not a subclass Y_l__m
%   Y = Y_lm(L, M) computes the spherical harmonic of orbital angular momentum
%   L and magnetic angular momentum M.
%
%   Y = Y_lm(__, Name, Value) specifies additional options:
%   - 'common': 'real' or 'complex'. Whether to compute the complex
%     spherical harmonics or their real part. Default is 'complex'.
%   - 'norm': Whether to normalize the result. Default is true.
%   - 'phase': Whether to include the Condon-Shortley phase. Default is true.

%   See also LEGENDRE,Y_L__M.

    properties(Access = private)
        m_real; % Stores the real magnetic quantum number (m)
        An = 0;
    end
    properties(Dependent)
        
    end
    methods
        function Y_lmObj = Y_lm(l,m,coe,n,options,options2)
            % Y_lm Constructor to initialize the spherical harmonic object.
            % Arguments:
            % l      - Azimuthal quantum number
            % m      - Magnetic quantum number
            % coe    - Coefficients for scaling
            % n      - Radius parameter (optional)
            % options - Additional options as name-value pairs
            arguments
                l   = 1; % Default Azimuthal quantum number
                m   = 0; % Default Magnetic quantum number
                coe = 1; % Default coefficient
                n   = 0; % Default radius
                options.common = true; % Default 
                options2.An = 0;
            end
            %
            PropertyCell = namedargs2cell(options);
            % See: https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
            % $$
            % \begin{cases}\frac{i}{\sqrt{2}}\left(Y_{\ell}^{-|m|}-(-1)^{m} Y_{\ell}^{|m|}\right) & \text { if } m<0 \\ Y_{\ell}^{0} & \text { if } m=0 \\ \frac{1}{\sqrt{2}}\left(Y_{\ell}^{-|m|}+(-1)^{m} Y_{\ell}^{|m|}\right) & \text { if } m>0\end{cases}
            % $$
            % Calculate the real coefficients for spherical harmonics
            absm = abs(m);
            if m < 0
                coeL(1)  =   +1i/sqrt(2);
                coeL(2)  =   -1i/sqrt(2)*(-1)^m;
            elseif m == 0
                coeL(1)  =  1/2;
                coeL(2)  =  1/2;
            elseif m > 0
                coeL(1)  =  +1/sqrt(2);
                coeL(2)  =  +1/sqrt(2)*(-1)^m;
            end
            % Adjust the input arguments if the spherical harmonics are inherited from another object
            if isa(l,'Y_l__m')
                mL = l.m;
                coeL = l.coe;
                l = l.l;
            else
                mL(1) = -absm;
                mL(2) =  absm;
                coeL = coeL * coe;
            end
            Y_lmObj = Y_lmObj@Y_l__m(l,mL,coeL,n,PropertyCell{:});
            for i = 1:length(Y_lmObj)
                Y_lmObj(i).m_real = m;
                Y_lmObj(i).An = options2.An;
            end
            
        end
    end
    methods % reload
         Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
    
         function SingleSum = InnerProduct_row(A_row, B_row, options)
            % Calculate the inner product for a pair of rows from HollowKnight objects A and B.
            % Options can specify symmetry, strictness, union, and squareness.
            arguments
                A_row
                B_row
                options.sym = true;        % Whether to use symbolic computation
                options.strict = false;    % Whether to enforce strict matching
                options.union = false;     % Whether to perform union operation
                options.square = true;     % Whether to enforce square matrix
            end
            % Initialize the sum with appropriate data type
            if options.sym
                SingleSum = sym(0);
            else
                SingleSum = 0;
            end
            if isrow(A_row) && isrow(B_row)
                if A_row(1).An == B_row(1).An
                    Tesseral_A = A_row.Tesseral;
                    Tesseral_B = B_row.Tesseral;
                    for i = 1:length(Tesseral_A{2})
                        iA_row = Tesseral_A{1}(i,:);
                        for j = 1:length(Tesseral_B{2})
                            jB_row = Tesseral_B{1}(j,:);
                            if all(iA_row == jB_row)
                                % Accumulate the product of coefficients for matching elements
                                SingleSum = SingleSum + Tesseral_A{2}(i) * Tesseral_B{2}(j);
                            end
                        end
                    end
                end
            end
        end
    end
    methods % disp
        function disp(YlmObj,options)
            % disp - Display method for the spherical harmonic object
            % Arguments:
            % options - Options for formatting the display (vpa, explicit, etc.)

            arguments
                YlmObj Y_lm;
                options.vpa = true; % Default to use variable precision arithmetic
                options.explicit = true; % Default to explicit formula
                options.cart = true; % Default to Cartesian format
            end
            optionsCell = namedargs2cell(options);
            if isscalar(YlmObj)
                disp(YlmObj.explicitformula(optionsCell{:}));
            else
                disp(YlmObj.formula(optionsCell{:}));
            end
        end
    end
    methods(Static)
    end
    methods %get
    end
    methods %math
        function Y_lmObj = SplusOper(Y_lmobj)
            Y_lmObj = Y_lmobj;
            for i = 1:numel(Y_lmobj)
                Y_lmObj(i).coe = Y_lmObj(i).coe * sqrt(Y_lmobj(i).s2_bar - Y_lmobj(i).Jz * (Y_lmobj(i).Jz + 1));
                Y_lmObj(i).Jz = Y_lmobj(i).Jz + 1;
            end
            Y_lmObj = Y_lmObj.contract();
        end
    end
end

