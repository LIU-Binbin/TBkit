function Bc = BerryCurvature_definition(A_1,A_2,para1,para2)
            Bc = diff(A_2,para1)-diff(A_1,para2);
        end