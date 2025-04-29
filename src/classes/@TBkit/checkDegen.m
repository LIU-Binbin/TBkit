function DegenPair = checkDegen(EIGEN,Accuracy)
%CHECKDEGEN Identify degenerate eigenvalue pairs
%
%   Syntax:
%       DegenPair = checkDegen(EIGEN,Accuracy)
%
%   Description:
%       Finds groups of degenerate eigenvalues within specified tolerance.
%
%   Inputs:
%       EIGEN    - Array of eigenvalues
%       Accuracy - Degeneracy tolerance (default=1e-8)
%
%   Output:
%       DegenPair - N×2 array of [start,end] indices for degenerate groups
arguments
    EIGEN {mustBeVector}
    Accuracy = 1e-8;
end
Nband= length(EIGEN);
DEGEN_init = (EIGEN(2:Nband)- EIGEN(1:(Nband-1)))<Accuracy;
tail = false;
DegenPair= [];
for i = 1:numel(DEGEN_init)
    if DEGEN_init(i)==1
        if  ~tail
            DEGEN_Pairlittle = i-1+1;
            tail = true;
        end
    elseif DEGEN_init(i)==0
        if tail
            DEGEN_Pairlittle = [DEGEN_Pairlittle,i-1+1];
            tail = false;
            DegenPair= [DegenPair;DEGEN_Pairlittle];
        end
    end
end
if tail
    DEGEN_Pairlittle = [DEGEN_Pairlittle,i+1];
    tail = false;
    DegenPair= [DegenPair;DEGEN_Pairlittle];
end
if ~tail
    return;
else
    error('?');
end
%DEGEN_char = reshape(char(string(DEGEN_init)),1,[]);
end