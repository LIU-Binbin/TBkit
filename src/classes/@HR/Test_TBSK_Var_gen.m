function EQ_list = Test_TBSK_Var_gen(testmode)
%TEST_TBSK_VAR_GEN Test function for TBSK_Var_gen
%
%   EQ_LIST = TEST_TBSK_VAR_GEN(TESTMODE) verifies the correctness of TBSK_Var_gen
%   by comparing generated symbolic expressions with expected forms.
%
%   Input:
%       testmode - Testing mode selector (default = 1)
%                  1: Tests with L_2 >= L_1
%                  Other: Tests all combinations
%
%   Output:
%       EQ_list  - List of logical equations showing equality comparisons
%
%   See also TBSK_VAR_GEN, HR.TBSK_Var_gen_single
arguments
    testmode =1 ;
end
if testmode == 1
    count = 0;
    for L_1 = 0:2
        for L_2 = L_1:2
            for m_1 = -L_1:L_1
                for m_2 = -L_2:L_2
                    count = count +1;
                    orb_sym1 = Ymlsym(L_1,m_1);
                    orb_sym2 = Ymlsym(L_2,m_2);
                    if orb_sym1 == sym(1)
                        orb_sym1 = 's';
                    end
                    if orb_sym2 == sym(1)
                        orb_sym2 = 's';
                    end
                    SymvarL = str2sym(...
                        strrep(...
                        strrep(...
                        strrep(['E_',char(orb_sym1),'_',char(orb_sym2)],'^',''),...
                        '*',''),...
                        ' - ','m')...
                        );

                    SymVarR = HR.TBSK_Var_gen_single(L_1,L_2,m_1,m_2,'sym_mode',1);
                    SymVarR_exchange = HR.TBSK_Var_gen_single(L_2,L_1,m_2,m_1,'sym_mode',1);
                    exchange_test = sym('Exchange') == simplify(SymVarR_exchange - SymVarR);
                    EQ(count,:) = [(SymvarL == SymVarR),exchange_test];
                end
            end
        end
    end
else
    count = 0;
    for L_1 = 0:2
        for L_2 = 0:2
            for m_1 = -L_1:L_1
                for m_2 = -L_2:L_2
                    count = count +1;
                    orb_sym1 = Ymlsym(L_1,m_1);
                    orb_sym2 = Ymlsym(L_2,m_2);
                    if orb_sym1 == sym(1)
                        orb_sym1 = 's';
                    end
                    if orb_sym2 == sym(1)
                        orb_sym2 = 's';
                    end
                    SymvarL = str2sym(...
                        strrep(...
                        strrep(...
                        strrep(['E_',char(orb_sym1),'_',char(orb_sym2)],'^',''),...
                        '*',''),...
                        ' - ','m')...
                        );

                    SymVarR = HR.TBSK_Var_gen_single(L_1,L_2,m_1,m_2,'sym_mode',1);
                    EQ(count,:) = (SymvarL == SymVarR);
                end
            end
        end
    end
end
EQ_list = EQ;
end
