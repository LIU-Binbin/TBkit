classdef lm < Spin
    % latter we may make Spin and lm jm as the subgroup of SO(3) or SU(2)
    % by Lie agbrea
    %  We choose the component lz and denote the common
    %  eigenstate of the operators l^2
    % \begin{aligned}
    % &L_{z}|l, m\rangle=m \hbar|l, m\rangle \\
    % &L^{2}|l, m\rangle=l(l+1) \hbar^{2}|l, m\rangle
    % \end{aligned}

    properties(Dependent,Access=private)
        l;
        m;
    end
    methods
        function lm_basis = lm(l,m,coe,options)
            % lm_basis = lm(1,0)
            %
            arguments
                l {lm.mustBeIntegerOrSpin(l)} = 1;
                m = [];
                coe = 1;
                options.orientation = [0 0 1];
            end
            optionsCell = namedargs2cell(options);
            if isa(l,'Spin')
                Spin_list= l;
                SIZE_Spin_list = size(Spin_list);
                l = reshape([Spin_list.J],SIZE_Spin_list);
                m = reshape([Spin_list.Jz],SIZE_Spin_list);
                coe = reshape([Spin_list.coe],SIZE_Spin_list);
            end
            if isempty(m)
                lz = l:-1:-l;
            else
                lz = m;
            end
            lm_basis = lm_basis@Spin(l,lz,coe,optionsCell{:});
            %lm_basis = lm_basis@Spin(l(1),lz(1),coe(1),optionsCell{:});
%             if max([length(Sz),length(l),length(coe)]) > 1
%                 lm_basis = repmat(lm_basis,[]);
%             end
%             for i = 2:max([length(lz),length(l),length(coe)])
%                 if length(lz) >1
%                     mi = i;
%                 else
%                     mi = 1;
%                 end
%                 if length(l) >1
%                     li = i;
%                 else
%                     li = 1;
%                 end
%                 if length(coe) >1
%                     ci = i;
%                 else
%                     ci = 1;
%                 end
%                 lm_basis(i) = lm(l(li),lz(mi),coe(ci),optionsCell{:});
%             end
        end
    end
    methods %get
        function l = get.l(lm_basis)
            l = lm_basis.J;
        end
        function m = get.m(lm_basis)
            m = lm_basis.Jz;
        end
    end
    methods
        function Lxm = Lx(lm_basis,options)
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lxm = Lx@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lxm = U*Lxm*U'/1i;
            else
                Lxm = Lx@Spin(lm_basis,optionsCell{:});
            end
        end
        function Lym = Ly(lm_basis,options)
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lym = Ly@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lym = U*Lym*U'/1i;
            else
                Lym = Ly@Spin(lm_basis,optionsCell{:});
            end
        end
        function Lzm = Lz(lm_basis,options)
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lzm = Lz@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lzm = U*Lzm*U'/1i;
            else
                Lzm = Lz@Spin(lm_basis,optionsCell{:});
            end
        end
    end
    methods(Static)
        function mustBeIntegerOrSpin(l)
            if isa(l,'Spin')
            else
                mustBeInteger(l);
            end
        end
    end
end
classdef lm < Spin
    %{
    This class represents the quantum states |l, m> where l is the 
    quantum number associated with L^2 and m is the quantum number 
    associated with Lz. 
    
    latter we may make Spin and lm jm as the subgroup of SO(3) or SU(2)
    by Lie agbrea
    
    The eigenvalues of these operators are:
    \begin{aligned}
    &L_{z}|l, m\rangle=m \hbar|l, m\rangle \\
    &L^{2}|l, m\rangle=l(l+1) \hbar^{2}|l, m\rangle
    \end{aligned}
    %}
    
    properties(Dependent,Access=private)
        l;
        m;
    end
    
    methods
        % Constructor for the lm class
        function lm_basis = lm(l,m,coe,options)
            %{
            This constructor creates an instance of the 'lm' class, which 
            represents a quantum state |l, m> where l is the quantum number 
            associated with L^2 and m is the quantum number associated with Lz. 
            The state is initialized with the given values of l, m, and coe 
            (coefficient).
            
            Inputs:
                l: Quantum number associated with L^2 (default 1).
                m: Quantum number associated with Lz (default is from -l to l).
                coe: Coefficient of the state (default is 1).
                options.orientation: Orientation for the spin state (default [0, 0, 1]).
            
            Outputs:
                lm_basis: An instance of the 'lm' class with specified properties.
            %}
            arguments
                l {lm.mustBeIntegerOrSpin(l)} = 1;
                m = [];
                coe = 1;
                options.orientation = [0 0 1];
            end
            optionsCell = namedargs2cell(options);
            if isa(l,'Spin')
                Spin_list= l;
                SIZE_Spin_list = size(Spin_list);
                l = reshape([Spin_list.J],SIZE_Spin_list);
                m = reshape([Spin_list.Jz],SIZE_Spin_list);
                coe = reshape([Spin_list.coe],SIZE_Spin_list);
            end
            if isempty(m)
                lz = l:-1:-l;
            else
                lz = m;
            end
            lm_basis = lm_basis@Spin(l,lz,coe,optionsCell{:});
        end
    end
    
    methods %get
        % Get method for l
        function l = get.l(lm_basis)
            %{
            This method returns the quantum number l associated with the
            Spin object.

            Inputs:
                lm_basis: An instance of the 'lm' class.

            Outputs:
                l: The quantum number l for the Spin object.
            %}
            l = lm_basis.J;
        end
        
        % Get method for m
        function m = get.m(lm_basis)
            %{
            This method returns the quantum number m associated with the
            Spin object.

            Inputs:
                lm_basis: An instance of the 'lm' class.

            Outputs:
                m: The quantum number m for the Spin object.
            %}
            m = lm_basis.Jz;
        end
    end
    
    methods
        % Method for Lx operator
        function Lxm = Lx(lm_basis,options)
            %{
            This method calculates the Lx operator for the given 'lm' object.
            It uses a transformation for the case when l = 1 and applies 
            the operator for other cases.

            Inputs:
                lm_basis: An instance of the 'lm' class.
                options: Options for the calculation (e.g., full or symmetric).

            Outputs:
                Lxm: The calculated Lx operator.
            %}
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lxm = Lx@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lxm = U*Lxm*U'/1i;
            else
                Lxm = Lx@Spin(lm_basis,optionsCell{:});
            end
        end
        
        % Method for Ly operator
        function Lym = Ly(lm_basis,options)
            %{
            This method calculates the Ly operator for the given 'lm' object.
            It uses a transformation for the case when l = 1 and applies 
            the operator for other cases.

            Inputs:
                lm_basis: An instance of the 'lm' class.
                options: Options for the calculation (e.g., full or symmetric).

            Outputs:
                Lym: The calculated Ly operator.
            %}
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lym = Ly@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lym = U*Lym*U'/1i;
            else
                Lym = Ly@Spin(lm_basis,optionsCell{:});
            end
        end
        
        % Method for Lz operator
        function Lzm = Lz(lm_basis,options)
            %{
            This method calculates the Lz operator for the given 'lm' object.
            It uses a transformation for the case when l = 1 and applies 
            the operator for other cases.

            Inputs:
                lm_basis: An instance of the 'lm' class.
                options: Options for the calculation (e.g., full or symmetric).

            Outputs:
                Lzm: The calculated Lz operator.
            %}
            arguments
                lm_basis lm;
                options.full = true;
                options.sym = true;
            end
            optionsCell = namedargs2cell(options);
            if lm_basis(1).l == 1 && options.full 
                Lzm = Lz@Spin(lm_basis,optionsCell{:});
                U = 1/sqrt(2) * [-1 0 1;-1i 0 -1i;0 sqrt(2) 0];
                Lzm = U*Lzm*U'/1i;
            else
                Lzm = Lz@Spin(lm_basis,optionsCell{:});
            end
        end
    end
    
    methods(Static)
        % Static method to ensure l is either an integer or a Spin object
        function mustBeIntegerOrSpin(l)
            %{
            This method checks if the input l is either an integer or a Spin 
            object.

            Inputs:
                l: The value to check.

            Outputs:
                None. Throws an error if the condition is not met.
            %}
            if isa(l,'Spin')
            else
                mustBeInteger(l);
            end
        end
    end
end
