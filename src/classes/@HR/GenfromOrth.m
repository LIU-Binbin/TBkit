function H_hr = GenfromOrth(H_hr,seed_r,seed_i,options)
% GENFROMORTH Generate Hamiltonian from orthogonal vector representation
%
%   H_hr = GENFROMORTH(H_hr,seed_r,seed_i,Accuracy,options) constructs
%   Hamiltonian coefficients from orthogonal vector representation.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object with vector hopping data
%       seed_r - Base name for real symbolic variables
%       seed_i - Base name for imaginary symbolic variables
%       Accuracy - Numerical accuracy for rounding
%       options - Structure with parameters:
%           fromCvectorL: Use combined vector (logical)
%
%   OUTPUT ARGUMENTS:
%       H_hr - Modified HR object with generated coefficients
%
%   NOTES:
%       - Converts vector hopping representation to symbolic coefficients
%       - Uses vpa for numerical accuracy control
%
%   SEE ALSO:
%       HR, vpa
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
    H_hr HR;
    seed_r = 'gamma__r_';
    seed_i = 'gamma__i_';
    options.Accuracy = 1e-6;
    options.fromCvectorL = true;
    
end
Accuracy = options.Accuracy;
nAccuracy = floor(-log(Accuracy)/log(10));
switch H_hr.Type
    case 'list'
        if H_hr.vectorhopping
            if  options.fromCvectorL
                CL= vpa(H_hr.CvectorL(:,1:rank(H_hr.CvectorL)),nAccuracy);
                SymVar_r = sym(seed_r,[size(CL,2),1],'real');
                CL = polish(CL,Accuracy);
                H_hr.HcoeL = CL(1:end/2,:)*SymVar_r + 1i*CL(end/2+1:end,:)*SymVar_r;
                H_hr.vectorhopping = false;
                return;
            end
            AL = vpa(H_hr.AvectorL,nAccuracy);
            BL = vpa(H_hr.BvectorL,nAccuracy);
            SymVar_r = sym(seed_r,[size(AL,2),1],'real');
            SymVar_i = sym(seed_i,[size(BL,2),1],'real');
            H_hr.HcoeL = AL*SymVar_r + 1i*BL*SymVar_i;
            H_hr.vectorhopping = false;
        end
end
end

function CL = polish(CL,Accuracy)
checkNum = [1,2,3,-1,-2,-3,1/2,-1/2,1/3,-1/3];
for IcheckNum = checkNum
    CL(abs(CL - IcheckNum) < Accuracy) = IcheckNum;
end

end
