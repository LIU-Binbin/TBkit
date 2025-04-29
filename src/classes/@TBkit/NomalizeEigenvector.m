function Eigenvetor=NomalizeEigenvector(Eigenvetor)
            Eigenvetor = ...
                Eigenvetor/((Eigenvetor'*Eigenvetor)^(1/length(Eigenvetor)));
        end