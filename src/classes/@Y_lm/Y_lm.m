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
    end
    properties(Dependent)
        
    end
    methods
        function Y_lmObj = Y_lm(l,m,coe,n,options)
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
            end
        end
    end
    methods % reload
        function Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
            % rotateinner - Rotates the spherical harmonics using Wigner D-matrix
            % Arguments:
            % obj - Current object
            % abc - Rotation angles [alpha, beta, gamma]
            % RightorLeft - Flag for right or left rotation
            % immproper, conjugate, antisymmetry - Optional parameters for symmetry handling

            arguments
                A
                abc
                RightorLeft
                immproper = true;
                conjugate = false;
                antisymmetry = false;
            end
            alpha = (abc(1));
            beta = (abc(2));
            gamma = (abc(3));
            A_L = Y_l__m(A.l);
            Ak = A.HollowMe;
            for i =  1:length(A_L)
                Ai = Y_lm(A_L(i));
                if A.l == 0
                    Ai.coe = 1;
                else
                    m1 = A.m;
                    m2 = Ai.m;
                    WignerD_single_element = (Y_l__m.d(A.l,m1,m2,beta));
                    Ai.coe = Ai.coe*exp(1i*RightorLeft*m1*alpha)*WignerD_single_element*exp(1i*RightorLeft*m2*gamma);
                    % for Y_l__m
                    Ai.coe = conj(Ai.coe);
                end
                if immproper
                    Ai.coe = (-1)^(A.l)*Ai.coe;
                end
                if ~conjugate
                    Ai.coe = Ai.coe * A.coe;
                else
                    Ai.coe = conj(Ai.coe) * A.coe; % check why the later cant use conj!
                end
                if Ai.coe ~= zeros(1,1,class(Ai.coe)) && ~isnan(Ai.coe)
                    Ak = [Ak,Ai];
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
        
    end
end

