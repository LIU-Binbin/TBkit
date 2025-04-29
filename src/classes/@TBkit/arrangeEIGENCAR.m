function EIGENCAR = arrangeEIGENCAR(EIGENCAR,REFCAR,method,opt)
%ARRANGEEIGENCAR Arrange eigenvalue array for band structure visualization
%
%   Syntax:
%       EIGENCAR = arrangeEIGENCAR(EIGENCAR,REFCAR,method,opt)
%
%   Description:
%       Rearranges eigenvalues to ensure smooth band connections by
%       matching nearest neighbors between k-points. Supports both
%       Hermitian and non-Hermitian systems.
%
%   Inputs:
%       EIGENCAR - Eigenvalue array (bands × k-points)
%       REFCAR   - Reference array for sorting (default=EIGENCAR)
%       method   - Sorting method:
%                  'NonHermitian' - For complex eigenvalues
%                  'PBAND_single' - For projected band structures
%       opt      - Options structure with fields:
%                  SelectL - Selected projections (PBAND_single only)
%                  Ecut    - Energy cutoff range [min,max]
%                  Disp    - Display debug info (logical)
%
%   Output:
%       EIGENCAR - Rearranged eigenvalue array
arguments
    EIGENCAR ;
    REFCAR = EIGENCAR;
    method {mustBeMember(method,{'NonHermitian','PBAND_single'})} = 'NonHermitian';
    opt.SelectL = [];
    opt.Ecut = [-3,3];
    opt.Disp = false;
end
switch  method
    case 'NonHermitian'
        %XList = real(REFCAR);
        %YList = imag(REFCAR);
        Nbands = size(EIGENCAR,1);
        Nk = size(EIGENCAR,2);
        EIGENCAR_OUT = EIGENCAR;
        NormelSeq = (1:Nbands).';
        % dsearchn
        for i = 2:Nk
            XListTemp = real(EIGENCAR_OUT(:,i-1:i));
            YListTemp = imag(EIGENCAR_OUT(:,i-1:i));

            PQ = [XListTemp(:,1),YListTemp(:,1)];
            P  = [XListTemp(:,2),YListTemp(:,2)];
            [k,dist] = dsearchn(P,PQ);

            if ~isequal(k,NormelSeq)
                if opt.Disp
                    disp(i);
                    disp(k);
                end
                for j = i:Nk
                    EIGENCAR_OUT(:,j) = EIGENCAR_OUT(k,j);
                end
            else

            end
        end
        EIGENCAR = EIGENCAR_OUT;
    case 'PBAND_single'
        %WEIGHTCAR  = REFCAR;
        Nbands = size(EIGENCAR,1);
        Nk = size(EIGENCAR,2);
        EIGENCAR_OUT = EIGENCAR;
        NormelSeq = (1:Nbands).';
        [NBANDS,NK,NPROJECTION]= size(REFCAR);
        if NBANDS ~= Nbands || Nk ~= NK
            error('REFCAR and EIGENCAR mismatch!');
        end
        if isempty(opt.SelectL)
            SelectL = 1:NPROJECTION;
        else
            SelectL = opt.SelectL;
        end
        Ecut =  opt.Ecut;
        % dsearchn
        for i = 2:Nk
            %XListTemp = real(EIGENCAR_OUT(:,i-1:i));
            %YListTemp = imag(EIGENCAR_OUT(:,i-1:i));
            PQ = reshape(REFCAR(:,i-1,:),NBANDS,NPROJECTION);
            P  = reshape(REFCAR(:,i,:),NBANDS,NPROJECTION);
            Eselect1 = Ecut(1) < EIGENCAR_OUT(:,i-1) & EIGENCAR_OUT(:,i-1) < Ecut(2) ;
            Eselect2 = Ecut(1) < EIGENCAR_OUT(:,i) & EIGENCAR_OUT(:,i) < Ecut(2) ;
            Eselect = find(Eselect1 & Eselect2);
            PQ = PQ(Eselect,SelectL);
            P = P(Eselect,SelectL);
            [k,~] = dsearchn(P,PQ);
            realk = [(1:min(Eselect)-1),k'+min(Eselect)-1,(max(Eselect)+1):NBANDS];
            if ~isequal(realk,NormelSeq)
                %disp(i);
                %disp(k);
                %for j = i:Nk
                %    EIGENCAR_OUT(:,j) = EIGENCAR_OUT(k,j);
                %end
                EIGENCAR_OUT(:,i:Nk) = EIGENCAR_OUT(realk,i:Nk);
            else

            end
        end
        EIGENCAR = EIGENCAR_OUT;
  end
end