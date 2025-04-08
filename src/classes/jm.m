classdef jm <Spin
    %{
    This class represents the quantum states denoted by |j, m>.
    These states are eigenstates of the operators J^2 and Jz. 
    The eigenvalues of these operators are:
    \begin{aligned}
    \boldsymbol{J}^{2}|j, m\rangle &=j(j+1) \hbar^{2}|j, m\rangle, j=0, \frac{1}{2}, 1, \frac{3}{2}, \ldots \\
    J_{z}|j, m\rangle &=m \hbar|j, m\rangle, m=-j,-j+1, \ldots, j-1, j
    \end{aligned}
    %}
    
    properties
        j
        m
    end
    
    methods
        % Constructor for the jm class
        function jm_basis = jm(j,jz,coe,options)
            %{
            This constructor creates an instance of the 'jm' class, 
            which represents a quantum state |j, m> where j is the 
            quantum number associated with J^2 and m is the quantum number 
            associated with Jz. The state is initialized with the given 
            values of j, jz, and coe (coefficient). 

            Inputs:
                j: Quantum number associated with J^2 (default 3/2).
                jz: Quantum number associated with Jz (default is from -j to j).
                coe: Coefficient of the state (default is 1).
                options.orientation: Orientation for the spin state (default [0, 0, 1]).
            
            Outputs:
                jm_basis: An instance of the 'jm' class with specified properties.
            %}
            arguments
                j = 3/2;
                jz = [];
                coe = 1;
                options.orientation = [0 0 1];
            end
            optionsCell = namedargs2cell(options);
            if isa(j,'Spin')
                Spin_list= j;
                j = [Spin_list.J];jz = [Spin_list.Jz];coe = [Spin_list.coe];
            end
            if isempty(jz)
                jz = j:-1:-j;
            else
                %jz = jz;
            end
            jm_basis = jm_basis@Spin(j,jz,coe,optionsCell{:});
        end
    end
    
    methods %get
        % Get method for j
        function j = get.j(lm_basis)
            %{
            This method returns the quantum number j associated with the
            Spin object.

            Inputs:
                lm_basis: An instance of the 'jm' class.

            Outputs:
                j: The quantum number j for the Spin object.
            %}
            j = lm_basis.J;
        end
        
        % Get method for m
        function m = get.m(lm_basis)
            %{
            This method returns the quantum number m associated with the
            Spin object.

            Inputs:
                lm_basis: An instance of the 'jm' class.

            Outputs:
                m: The quantum number m for the Spin object.
            %}
            m = lm_basis.Jz;
        end
    end
end
