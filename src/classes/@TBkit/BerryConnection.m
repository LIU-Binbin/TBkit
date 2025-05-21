        function F = BerryConnection(W1,W2)
            %             Norb1 = size(W1,2);
            %             Norb2 = size(W2,2);
            %             for n = 1:Norb1
            %                 for m = 1:Norb2
            %                     F(n,m) = W1(:,n)'* W2(:,m);
            %                 end
            %             end
            F = W1'*W2;
        end