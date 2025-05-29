function [FITobj,SubsIndexL] = fitprepare(FITobj,dH,SelectVectorL)
arguments
    FITobj
    dH = []
    SelectVectorL = -1;
end
if SelectVectorL == -1
    try
        SelectVectorL = FITobj.Line_000;
    catch

    end
end
if isempty(dH)
    FITobj.HnumL = zeros(size(FITobj.HcoeL));
    num_seqL = ones(1,numel(FITobj.HnumL ));
    SubsIndexL = num_seqL;
    ntotal = numel(FITobj.HnumL);
    for i = 1:numel(FITobj.HnumL)
        fprintf('%d/%d\n',i,ntotal);
        try
            FITobj.HnumL(i) = double(FITobj.HcoeL(i));
            % num_seqL(i) = 0;
        catch
            SubsIndexL(i) = 0;
        end

    end
    SubsIndexL = find(SubsIndexL);
    % FITobj.HnumL(num_seqL) = FITobj.HcoeL(num_seqL);
    FITobj.HcoeL  = FITobj.HcoeL(SubsIndexL);
else
    SelectM = find(logical(dH ~= sym(0)));
    size_dH = size(dH);
     % = ind2sub(size_dH,SelectM);
    [rowL,colL,page] = ind2sub(size_dH,SelectM);
    pageL = SelectVectorL(page);
    SubsIndexL = sub2ind(size(FITobj.HnumL),rowL,colL,pageL);
    FITobj.HcoeL =dH(SelectM);
    
end
end