classdef group
    %% Group operation of discontinuous groups.
    % This class implements operations on groups including identity element, 
    % inverse elements, closure, abelian property, and more.

    %% Description of the class:
    %     U : double (optional)
    %         The unitary action on the Hilbert space.
    %         May be None, 
    %
    %     * Late 1700s- Joseph-Louis Lagrange (1736-1813) 利用置换的概念，理解了三次和四次方程为什么有解。（Paolo Ruffini利用同样的思想证明了一般的五次或以上方程没有根式解）
    %     * Early 1800s-Évariste Galois (killed in a duel in 1832 at age 20), and Niels
    %     * Abel (died in 1829 at age 26 of TB) 阐述了代数方程的可解性和置换群的联系，真正利用群的概念解决了这个难题。
    %     * 第一个正式的群的概念由Cayley于1854提出。Cayley, Hamilton, and Sylvester 引进了矩阵的表达。1890年Fedorov发展了各种对称性操作的数学方法，证明了晶体有且仅有230种空间群，次年Schönflies也独立的证明了这个结果。1920年群论开始应用于物理、化学等领域。里面还有很多著名数学家和物理学家的贡献，不一一而举。
    properties
        U; % The unitary action on the Hilbert space.
    end

    properties(GetAccess = protected, Hidden = true)
        isIdentity = false; % Flag for identity element.
        isClosed = false;   % Closure flag.
    end

    methods
        % Constructor to initialize the group with a unitary matrix.
        function obj = group(U)
            % If no unitary matrix is provided, initialize as NaN.
            arguments
                U = nan;
            end
            obj.U = U;
        end
    end

    %% Group Defination
    % * 单位元：∃𝑒∈𝑆,∀𝑥∈𝑆,𝑒⋅𝑥=𝑥⋅𝑒=𝑥;
    % Identity element
    % There exists an element 𝑒⋅𝑥=𝑥⋅𝑒=𝑥; It is called the identity element of the group.
    methods(Static)
        function E = identity()
            % Return the identity element of the group.
            E = group();
            E.isIdentity = true;
        end
    end
    methods
        % E
        % Return the identity element of the same structure as the group.
        function E = E(group)
            if ~isnan(group.U)
                group.U = eye(length(group.U), class(group.U));
            end
            E = group;
            E.isIdentity = true;
        end
    end
    %  Associativity
    % For all a, b, c in G, one has (a・b)・c=a・(b・c). Alreay imply in ();
    % Inverse element
    % Close: ∀𝑥,𝑦∈𝑆,𝑥⋅𝑦∈𝑆;  imply in * ;
    %% Group operation overloads
    methods
        % ==
        function TrueOrFalse = eq(G1,G2)
            % Check whether two elements or groups are identity
            %
            m = length(G1);
            n = length(G2);
            if m == 1 && n ==1
                if G1.isIdentity && G2.isIdentity
                    TrueOrFalse = true;
                    return;
                end
                if isequal(G1.U,G2.U)
                    TrueOrFalse = true;
                else
                    TrueOrFalse = false;
                end
                return;
            elseif m == 0 && n == 0
                TrueOrFalse =true;
                return;
            elseif m == 0 || n == 0
                TrueOrFalse = false;
                return
            elseif m == 1
                TrueOrFalse(n) = false;
                for i = 1:n
                    TrueOrFalse(i) = eq(G1,G2(i));
                end
                return;
            elseif n == 1
                TrueOrFalse(m) = false;
                for i = 1:m
                    TrueOrFalse(i) = eq(G1(i),G2);
                end
                return;
            else
                % unique need sort( sort need overload lt(less than))
                % G1 = unique(G1);
                % G2 = unique(G2);
                %
                % m = length(G1);
                % n = length(G2);
                if m == n
                    TrueOrFalse = false(size(G1));
                    for i =1:m
                        TrueOrFalse(i) = eq(G1(i),G2(i));
                    end
                else
                    % generate_group ?
                    TrueOrFalse =false;
                end
            end
        end
        % ~=
        function result = ne(G1, G2)
            % Inequality check for two groups
            result = ~(G1 == G2);
        end
        % <
      function TrueOrFalse = lt(G1,G2)
            % Less than check based on unitary matrix comparison
            % Sort group elements:
            if G1.isIdentity && ~G2.isIdentity
                TrueOrFalse = true;
                return;
            elseif ~G1.isIdentity && G2.isIdentity
                TrueOrFalse = false;
                return;
            else
            end
            L1 = logical(real(G1.U(:)) < real(G2.U(:)));
            B1 = logical(~(real(G1.U(:)) == real(G2.U(:))));
            L2 = logical(imag(G1.U(:)) < imag(G2.U(:)));
            B2 = logical(~(imag(G1.U(:)) == imag(G2.U(:))));
            for i =1:length(B1)
                if B1(i)
                    TrueOrFalse = L1(i);
                    return;
                end
            end
            for i =1:length(B2)
                if B2(i)
                    TrueOrFalse = L2(i);
                    return;
                end
            end
            TrueOrFalse  = false;
        end
        % >
        % Greater than check based on unitary matrix comparison
        function TrueOrFalse = gt(G1,G2)
            TrueOrFalse = lt(G2,G1);
        end
        % inv
        % Inverse of the group
        function G = inv(G1)
            % Invert Operator
            % Standerd definition
            % need reload /
            if isscalar(G1)
                G =G1;
                if isnan(G.U)
                else
                    G.U = inv(G.U);
                end
            else
                for i = 1:length(G1)
                    G1(i) = inv(G1(i));
                end
            end
        end
        % uminus
        % Unary minus operation
        function G = uminus(G)
            for i = 1:numel(G)
                G(i).U = -G(i).U;
            end
        end
        % +
        % Add two groups together (combine unique elements)
        function G = plus(G1,G2)
            G = unique([G1,G2]);
        end
        % -
        % Subtract one group from another
        function G = minus(G1,G2)
            % unique need sort( sort need overload lt(less than))
            % Remove elements of G2 from G1 and return the result.
            G = setdiff(G1, G2, 'sorted');
        end
        % /
        function G = mrdivide(G1, G2)
            % Perform the group division operation (G1 / G2).
            % First, ensure both groups have unique elements, then check the condition on lengths.

            % Ensure the elements of both groups are unique
            G1 = unique(G1);
            G2 = unique(G2);

            % Get the number of elements in each group
            m = length(G1);
            n = length(G2);

            % Initialize an empty result
            G = [];

            % Check if the division condition (m is a multiple of n) is satisfied
            if mod(m, n) == 0
                % Perform the quotient operation if the condition holds
                G = quotient(G2, G1);
            else
                % Return an empty result if division is not valid
                G = [];
            end
        end
        % ./
        function G3 = rdivide(G1, G2)
            % Perform division operation for two groups or matrices.
            % The division logic varies depending on whether G1 and G2 are identity elements or not.

            % Ensure both G1 and G2 have unique elements
            G1 = unique(G1);
            G2 = unique(G2);

            % Handle the case when both G1 and G2 are scalar-like or identity elements
            if isscalar(G1) && isscalar(G2)
                % Case when G1 is the identity element
                if G1.isIdentity
                    G3 = G2;
                    G3.U = inv(G2); % Inverse of G2
                    return;
                    % Case when G2 is the identity element
                elseif G2.isIdentity
                    G3 = G1;
                    return;
                    % Case for general multiplication and inverse
                else
                    G3 = G1;
                    G3.U = G1 .* G2.inv(); % Inverse of G2, multiplied by G1
                end
            else
                % Error handling for unsupported operations (not elementary operator)
                error('.* or ./ is an elementary operator and can only be used with scalar or identity elements');
            end
        end
        % *
        function G3 = mtimes(G1, G2)
            % mtimes - Define a matrix multiplication function with optimized implementation
            % Input:
            %   G1 - First input group
            %   G2 - Second input group
            % Output:
            %   G3 - group containing unique elements from the product of two input group

            % 1. Remove duplicates and sort input groups
            % Use the unique function to remove duplicate elements and sort them in ascending order
            G1 = unique(G1);
            G2 = unique(G2);

            % 2. Get group lengths
            m = length(G1); % Number of elements in G1
            n = length(G2); % Number of elements in G2

            % 3. Initialize result array
            % Use repmat to create an array of the same size as the result, initialized with the first element of G1
            newgroup = repmat(G1(1), 1, m*n);

            % 4. Calculate all products
            count = 0; % Initialize counter
            for i = 1:m
                for j = 1:n
                    count = count + 1;
                    newgroup(count) = G1(i) .* G2(j);
                end
            end

            % 5. Remove duplicates to get the final result
            G3 = unique(newgroup);
        end
        % .*
        function G3 = times(G1, G2)
            % Perform multiplication for two groups or matrices.
            % The operation accounts for identity elements and simplifies symbolic results.

            % Case where both G1 and G2 are scalar-like or identity elements
            if isscalar(G1) && isscalar(G2)
                % Case when G1 is the identity element
                if G1.isIdentity
                    G3 = G2;
                    G3.U = inv(G2); % Inverse of G2
                    return;
                    % Case when G2 is the identity element
                elseif G2.isIdentity
                    G3 = G1;
                    return;
                    % General case for multiplication of the U components
                else
                    G3 = G1;
                    G3.U = G1.U * G2.U; % Multiply U components

                    % Simplify if the result involves symbolic elements
                    if isa(G3.U, 'sym')
                        G3.U = simplify(G3.U);
                    end
                end
            else
                % Error handling when the operation is not defined for these types
                error('.* is an elementary operator and can only be used with scalar or identity elements');
            end
        end
        % .^
        function G3 = power(G1, n)
            % Calculate the power of the group G1 raised to the integer exponent n.
            % The operation accounts for both positive and negative exponents.

            % Initialize the result with the identity element of G1
            G3 = G1.E(); % Assuming E() returns the identity element of G1

            % Loop to multiply G1 by itself |n| times
            for i = 1:abs(n)
                G3 = G3 .* G1; % Perform multiplication with G1
            end

            % If the exponent is negative, take the inverse of the result
            if n < 0
                G3 = G3.inv(); % Assuming inv() computes the inverse of G3
            end
        end
        % sort
        function [group_list, indSort] = sort(group_list)
            % Sort the elements of the group_list using insertion sort.
            % The function returns the sorted list and the indices of the sorted elements.

            % Call the insertsort method from the group object to perform sorting
            [group_list, indSort] = insertsort(group_list);
        end
    end
    %% Group theory
    methods
        function order = order(SymOper)
            % Returns the order (length) of the SymOper object.
            % This method calculates the length of the SymOper, which
            % represents the number of elements in the operator.

            order = length(SymOper);  % Calculate and return the length of the object
        end
    end

    methods
        % Closure check: Determines if the group is closed under the group operation.
        function TrurOrFalse = closure(group)
            TrurOrFalse = true;
            if all([group.isClosed])
                return
            else
                for i = 1:numel(group)
                    if find(ismember(group(i) * group, group) == 0)
                        TrurOrFalse = false;
                        return;
                    end
                end
            end
        end
        % Abelian check: Checks if the group is Abelian (commutative).
        function TrurOrFalse = abelian(group)
            TrurOrFalse = true;
            for Gi = group
                for Gj = group
                    if Gi*Gj ~=Gj*Gi
                        TrurOrFalse =false;
                        return;
                    end
                end
            end
        end
        % Generate a group from generators
        function group = generate_group(gens)
            % Parameters:
            % gens : iterable of group generators
            %
            % Returns:
            % group : generated group from the generators

            % need reload minus
            gens = unique(gens); % Remove duplicate generators
            % here we keep all the elements generated so far
            group = gens; % Initialize group with generators
            % these are the elements generated in the previous step
            oldgroup = gens; % Store old group elements
            %            fprintf('group muplicity: %d\n',group.order);
            while true
                % Generate new elements by multiplying with generators
                newgroup= unique(oldgroup * gens);
                %                fprintf('newgroup muplicity: %d\n',newgroup.order);
                % only keep those that are new
                newgroup = newgroup - group; % Remove previously generated elements 
                %                fprintf('newgroup muplicity (after -): %d\n',newgroup.order);
                % if there are any new, add them to group and set them as old
                if ~isempty(newgroup)
                    group = unique([group,newgroup]); % Add new elements to the group
                    oldgroup = newgroup; % Update old group
                    % if there were no new elements, we are done
                else
                    break; % If no new elements, stop
                end
            end
        end
        % Find generators of the group
        function generator = generator(group,options)
            arguments
                group
                options.fast = true; % Option to use fast(Do not check closure)
            end
            if ~options.fast
                % Generate group if not closed
                if ~group.closure
                    group = group.generate_group();
                end
            end
            % Get group order
            Order = group.order();
            pool = 1:Order;
            findit = false;
            count = 0;

            % Check for minimal generators
            for i = pool
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    jgroup = group(ChooseL(j,:));
                    if isequal(length(jgroup.generate_group),Order)
                        findit = true;
                        count = count+1;
                        generator{count} = jgroup;
                        ngenerator = i;
                    end
                end
                if findit
                    break; % Exit loop once generator is found
                end
            end
        end
        % Rank of the group: Find the minimal number of generators required
        function [Rank,jgroup] = rank(group)
            if ~group.closure
                group = group.generate_group();
            end
            pool = 1:numel(group);
            for i = pool
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    jgroup = group(ChooseL(j,:));
                    if isequal(jgroup.generate_group,group)
                        Rank = i;
                        return;
                    end
                end
            end
        end
        % Generate all subgroups of the group, including trivial and itself
        % subgroup
        function subgroup = subgroup(group)
            %     Generate all subgroups of group, including the trivial group
            %     and itself.
            %
            %     Parameters
            %     ----------
            %     group : set of PointGroupElement
            %         A closed group as set of its elements.
            %
            %     Returns
            %     -------
            %     subgroups : dict of forzenset: set
            %         frozesets are the subgroups, sets are a generator set of the
            %         subgroup.
            %
            %     # Make all the cyclic subgroups generated by one element
            % hard coding , solve it if needed
            count = 0;
            count = count + 1;
            subgroup{1} = group(1).E; % Trivial subgroup
            count = count + 1;
            subgroup{2} = group; % Whole group as a subgroup
            if isprime(group.order())
                % 𝑝(𝑝是素数)阶群𝐺均是Abel群,且均同构于整数模𝑝的加法群ℤ𝑝.
                %for igroup = (group(subgroup{1} ~= group))
                %    count = count + 1;
                %    subgroup{count} = igroup;
                %end
                % Prime-order groups are Abelian and isomorphic to Z_p
                return;
            else
                % Use Lagrange's Theorem to find possible subgroups
                % Lagrange定理 G的每个子群的阶数都是𝐺的阶数的因数.
                plist = SolvingFactor(group.order());
                plist(1:2) = [];
                pool = 1:numel(group);
                for i = plist
                    ChooseL = nchoosek(pool,i);
                    for j = 1:size(ChooseL,1)
                        igroup = group(ChooseL(j,:));
                        if igroup.closure
                            count = count + 1;
                            subgroup{count} = igroup;
                        end
                    end
                end
            end
        end
        % Normal subgroups: Find all normal subgroups of the group
        function NormalSubgroup = normalsubgroup(group)
            count = 0;
            count = count + 1;
            NormalSubgroup{count} = group(1).E; % Trivial subgroup
            count = count + 1;
            NormalSubgroup{count} = group; % Whole group
            if isprime(group.order()) % Handle prime-order groups
                % 𝑝(𝑝是素数)阶群𝐺均是Abel群,且均同构于整数模𝑝的加法群ℤ𝑝.
                return;
            end
            % 若G是交换群, 则G的所有子群都是正规子群。
            % If the group is Abelian, all subgroups are normal
            if group.abelian
                NormalSubgroup = subgroup(group);
            end
            % Further handling using conjugation and normal subgroup properties
            % Lagrange定理 + 正规子群定义
            ConjugationClassifyCollection  = conjugationseparation(group);
            nC = length(ConjugationClassifyCollection);
            pool = 1:nC;
            plist = SolvingFactor(group.order());
            plist(1:2) = [];
            for i = 1:nC
                ChooseL = nchoosek(pool,i);
                for j = 1:size(ChooseL,1)
                    ConjugationClassifyCollectionTmp = ConjugationClassifyCollection(ChooseL(j,:));
                    SubgroupTmp = [ConjugationClassifyCollectionTmp{:}];
                    nSi = length(SubgroupTmp);
                    if ismember(nSi,plist)
                        if SubgroupTmp.closure()
                            count = count + 1;
                            NormalSubgroup{count} = SubgroupTmp;
                        end
                    else
                        continue;
                    end
                end
            end
        end
        % Separate a group into left or right cosets of a subgroup
        function CosetClassify  = cosetseparation(H,G,leftorright)
            arguments
                H
                G
                leftorright {mustBeMember(leftorright,{'l','r'})} =  'l'; % Default is left coset
            end
            % Ensure H is a subgroup of G
            if ~all(ismember(H, G))
                error('H must be a subgroup of G.');
            end

            % Initialize the pool of elements not yet assigned to a coset
            pool = G-H;
            count = 1;
            CosetClassify{1}=H; % First coset is the subgroup itself
            while ~isempty(pool)
                % Select the first element from the pool
                iG = pool(1);
                pool = pool(2:end); % Remove the selected element from the pool
                switch leftorright
                    case 'l'
                        Coset = iG * H; % Compute left coset
                    case 'r'
                        Coset = H * iG; % Compute right coset
                    otherwise
                        error('Invalid coset type. Use ''l'' for left or ''r'' for right.');
                end
                % Add the coset to the classification
                CosetClassify{count} = Coset;

                % Remove the coset elements from the pool
                pool = pool - Coset;
            end
        end
        % is Normal subgroup？ 
        function TrueOrFalse = isNSG(H, G)
            % isNSG: Determine if subgroup H is normal in group G
            % Input:
            %   H: subgroup, represented as a struct with elements and closure property
            %   G: group, represented as a struct with group operation and inverse
            % Output:
            %   TrueOrFalse: logical indicating if H is normal in G

            % Assume H is normal initially
            TrueOrFalse = true;

            % Check if H is closed under the group operation
            if ~H.closure
                TrueOrFalse = false;
                return;
            end

            % Initialize pool with elements of H
            pool_origin = H;
            pool = pool_origin;

            % Process each element in the original pool
            for i = 1:length(pool_origin)
                Hi = pool_origin(i);

                % Skip if Hi has already been processed
                if ~ismember(Hi, pool)
                    continue;
                end

                % Compute the conjugation class of Hi in G
                ConjugationClass = conjugation(Hi, G);

                % Check if all conjugates are in H
                if ~all(ismember(ConjugationClass, H))
                    TrueOrFalse = false;
                    return;
                end

                % Remove the conjugation class from the pool
                pool = pool - ConjugationClass;
            end
        end
        % Quotient group: Check if H is a normal subgroup of G and form the quotient group
        function factor_group = quotient(H, G, leftorright)
            arguments
                H
                G
                leftorright = 'l';  % Default to left coset quotient group
            end

            if H.isNSG(G)  % Check if H is normal in G
                factor_group = cosetseparation(H, G, leftorright);
            else
                factor_group = [];  % If H is not normal, return empty
            end
        end
         % Conjugation: Check if G1 is conjugate to G2 in G
        function TrueOrFalse = conjugate(G1,G2,G)
            TrueOrFalse = false;
            for i = 1:numel(G)
                Gi = G(i);
                if G1 == Gi*G2*Gi.inv()
                    TrueOrFalse = true;
                    return;
                end
            end
        end
        % Conjugation class: Find conjugates of G1 in G
        function ConjugationClass  = conjugation(G1,G)
            ConjugationClass = G1;
            for i = 1:numel(G)
                Gi = G(i);
                if conjugate(G1,Gi,G)
                    ConjugationClass = [ConjugationClass,Gi];
                end
            end
            ConjugationClass = unique(ConjugationClass);
        end
        % Conjugation separation: Separate G into conjugacy classes
        function ConjugationClassify  = conjugationseparation(G)
            ConjugationClassify{1} = G;
            ConjugationClassify{1}(:) = [];
            count = 0;
            pool_origin = G;
            pool = pool_origin;
            for Gi = pool_origin
                if ismember(Gi,pool)
                    count = count +1 ;
                    ConjugationClassify{count} = conjugation(Gi,G);
                    pool = pool - ConjugationClassify{count};
                end
            end
        end
        % 
    end
    %% tools
    methods(Static)

    end
        % algorithm
       
        %
 
end