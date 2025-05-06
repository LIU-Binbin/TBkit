function Output = pretty(SpinObj,options)
arguments
    SpinObj Spin
    options.output char{mustBeMember(options.output,{'str','sym','latex'})} = 'str';
    options.formatspec =  '%5.3f';
end
if isempty(SpinObj)
    Output = "";
    return;
end
optionsCell = namedargs2cell(options);
rightbraket = '\x27E9';
if isscalar(SpinObj)
    if SpinObj.hollow
        Output = "";
        return;
    end
    switch options.output
        case 'str'
            if SpinObj.coe == 1 || SpinObj.coe == sym(1)
                strcoe = "";
            elseif SpinObj.coe == -1 || SpinObj.coe == sym(-1)
                strcoe = "-";
            elseif ~isreal(SpinObj.coe)
                strcoe = string(sym(SpinObj.coe));
            elseif isa(SpinObj.coe,'sym')
                strcoe = string(sym(SpinObj.coe));
            else
                strcoe = num2str((SpinObj.coe),options.formatspec);
            end
            Output = strcoe +"|"+string(sym(SpinObj.J))+","+string(sym(SpinObj.Jz))+rightbraket+" ";
        case 'latex'
        case 'sym'
    end
elseif length(SpinObj) > 1
    switch options.output
        case 'str'
            for i = 1:size(SpinObj,1)
                Output(i,1) = pretty(SpinObj(i,1),optionsCell{:});
                for j =2:size(SpinObj,2)
                    Output(i,1) = Output(i,1) +"+ " +  pretty(SpinObj(i,j),optionsCell{:});
                end
                Output(i,1) = Output(i,1)+"\n";
            end
        case 'latex'
        case 'sym'
    end
else
end
end
