function Tesseralexpansion = Tesseral(Y_lmObj,options)
arguments
    Y_lmObj;
    options.sym = true;
end
% take care! Y_lm is actually Y_l__m with a vectorized
% combination
% we use Tesseral get real SH coe
% Y_{\ell}^{m}= \begin{cases}\frac{1}{\sqrt{2}}\left(Y_{\ell|m|}-i Y_{\ell,-|m|}\right) & \text { if } m<0 \\ Y_{\ell 0} & \text { if } m=0 \\ \frac{(-1)^{m}}{\sqrt{2}}\left(Y_{\ell|m|}+i Y_{\ell,-|m|}\right) & \text { if } m>0\end{cases}
for i = 1:size(Y_lmObj,1)       % l m coe
    Y_lmObjtmp = Y_lmObj(i,:);
    comparerowL = [];
    if options.sym
        sumrowL = sym([]);
    else
        sumrowL = [];
    end
    for j = 1:length(Y_lmObjtmp)
        l = Y_lmObjtmp(j).l;
        m = Y_lmObjtmp(j).m;
        %
        if m < 0
            coe1 = 1/sqrt(2);
            coe2= -1i/sqrt(2);
        elseif m>0
            coe1 = 1/sqrt(2)*(-1)^m;
            coe2 = 1i/sqrt(2)*(-1)^m;
        else
            coe1 = 1/2;
            coe2 = 1/2;
        end
        comparerowL = [comparerowL;l abs(m);l -abs(m)];
        sumrowL = [sumrowL;coe1*Y_lmObjtmp(j).coe;coe2*Y_lmObjtmp(j).coe];
    end
    [Tesseralexpansion{i,1},Tesseralexpansion{i,2}] = Y_l__m.generalcontractrow(comparerowL,sumrowL);
end
end
