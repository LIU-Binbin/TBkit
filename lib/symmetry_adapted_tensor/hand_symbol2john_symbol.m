function john_symbol = hand_symbol2john_symbol(String)
% |      intrinsic AHC |         a{V2} | *RevModPhys. 82.1539*    |
% | -----------------: | ------------: | ------------------------ |
% |  intrinsic 2nd AHC |        a{V2}V | *PhysRevLett.127.277202* |
% |  intrinsic 3rd AHC |     a{V2}[V2] | *PhysRevB.107.075411*    |
% |      extrinsic AHC |          {V2} |                          |
% |  extrinsic 2nd AHC |         {V2}V | *PhysRevLett.115.216806* |
% |  extrinsic 3rd AHC |        {V2}V2 | *PhysRevB.105.045118*    |
% |                    |               |                          |
% |      intrinsic PHC |        e{V2}V | *PhysRevLett.132.056301* |
% |  intrinsic 2nd PHC |       e{V2}V2 | *PhysRevLett.130.126303* |
% |  extrinsic 2nd PHC |      ae{V2}V2 | *PhysRevB.108.075155*    |
% |                    |               |                          |
% |      intrinsic SHC |           eV3 | *RevModPhys. 87.1213*    |
% |  intrinsic 2nd SHC |           eV4 | *PhysRevLett.134.056301* |
% |      extrinsic SHC |          aeV3 |                          |
% |                    |               |                          |
% |               CISP |          aeV2 |                          |
% | intrinsic 2nd CISP |          aeV3 | *PhysRevLett.129.086602* |
% | extrinsic 2nd CISP |           eV3 | *PhysRevLett.130.166302* |
    switch String
        case {'AHC','1AHC','1AHC_i'}
            john_symbol = 'a{V2}';
        case {'2AHC','2AHC_i'}
            john_symbol = 'a{V2}V';
        case {'2AHC_e'}
            john_symbol = '{V2}V';
    end
end