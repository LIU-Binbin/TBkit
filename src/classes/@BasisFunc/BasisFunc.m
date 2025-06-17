classdef BasisFunc < HollowKnight
properties
BForb = [];
BFnum =[];
BFuncL = [];
spin = [];
end
properties(Hidden,Dependent)
hollow;
end
properties(Access = private,Dependent)
FuncNum ;
BFclassL ;
end
properties(Hidden)
parity;
Rm = [1 0 0;0 1 0;0 0 1];
end
%% construction inner
methods
    function BasisFunction = BasisFunc(BFuncLOrigin,spin,BFnum,coe,orb,options)
        % BasisFunction = BasisFunc({Spin(1/2)})
        % BasisFunction = BasisFunc(Spin(1/2))
        arguments
            BFuncLOrigin = Qnum().HollowMe;
            spin  Spin = Spin(0,0) ;
            BFnum = [];
            coe = sym([]);
            orb = [];
            options.raw = true
            options.fast = true;
        end
        optionsCell = namedargs2cell(options);
        if isempty(BFuncLOrigin)
            BasisFunction = BasisFunc.empty([1 0]);
            return;
        end
        if isa(BFuncLOrigin,'double')
            SizeBasisFunc = BFuncLOrigin;
            if isempty(BFuncLOrigin)
                BasisFunction = BasisFunc.empty([1 0]);
            elseif length(SizeBasisFunc) == 1
                BasisFunction = BasisFunc([SizeBasisFunc 1]);
                return;
            else
                BasisFunction1 = BasisFunc();
                BasisFunction = repmat(BasisFunction1,SizeBasisFunc);
                return;
            end
        elseif isa(BFuncLOrigin,'TBkit')
            BasisFunction = BasisFunc.BasisFunction(BFuncLOrigin);
            return;
        end
        if isa(BFuncLOrigin,'Spin') ||isa(BFuncLOrigin,'jm')
            spin = BFuncLOrigin;
            BFuncLOrigin = sym(1);
        end
        if isempty(BFnum)                % coe
            BFnum = 1;
        end
        if isempty(coe)
            coe = sym(1);
        end
        BasisFunction(1,1).BFuncL = BFuncLOrigin(1,1);
        BasisFunction(1,1).spin = spin(1,1);
        BasisFunction(1,1).BFnum = BFnum(1,1);
        BasisFunction(1,1).coe = coe(1,1);
        if ~isempty(orb)
            BasisFunction(1,1).BForb = orb(1,:);
        end
        if isscalar(BFuncLOrigin) && isscalar(spin)
            return;
        end
        sizeBFuncL = size(BFuncLOrigin);
        sizespin = size(spin);
        if isequal(sizeBFuncL,sizespin)
        elseif isvector(spin) && isvector(BFuncLOrigin)
            spin = repmat(spin,sizeBFuncL);
            BFuncLOrigin = repmat(BFuncLOrigin,sizespin);
        else
            try
                if sizeBFuncL(2) > sizespin(1)
                    spin = repmat(spin,sizeBFuncL./sizespin);
                else
                    BFuncLOrigin = repmat(BFuncLOrigin,sizespin./sizeBFuncL);
                end
            catch
                error('wrong input!')
            end
        end
        sizeBF = size(BFuncLOrigin);
        if ~isequal(sizeBF,size(BFnum))
            BFnum = repmat(BFnum,sizeBF./size(BFnum));
        end
        if ~isequal(sizeBF,size(coe))
            coe = repmat(coe,sizeBF./size(coe));
        end
        BasisFunction = repmat(BasisFunction(1,1),sizeBF);
        for i = 1:numel(BasisFunction)
            BasisFunction(i) = BasisFunc(BFuncLOrigin(i),spin(i),BFnum(i),coe(i),optionsCell{:});
        end
        if ~isempty(orb)
            norbL = size(orb,1);
            if norbL == 1
                for i = 1:numel(BasisFunction)
                    BasisFunction(i).BForb = orb(1,:);
                end
            elseif norbL == numel(BasisFunction)
                for i = 1:numel(BasisFunction)
                    BasisFunction(i).BForb = orb(i,:);
                end
            elseif norbL == size(BasisFunction,1)
                for i = 1:size(BasisFunction,1)
                    for j = 1:size(BasisFunction,1)
                        BasisFunction(i,j).BForb = orb(i,:);
                    end
                end

            else

            end
        end
    end
end
%% construction outer
methods(Static)
    function BasisFunction = BasisFunction(TBkitobj)
        switch class(TBkitobj)
            case {'TBkit','HR','Htrig','HK'}
                % Qnum
                [BFuncLOrigin,S,SzL] = Qnum.QnumL(TBkitobj);
                spinL = Spin(S,SzL);
                BasisFunction = BasisFunc(BFuncLOrigin,spinL,1,1,TBkitobj.orbL);
                for i = 1:numel(BasisFunction)
                    BasisFunction(i).Rm = TBkitobj.Rm;
                end
        end

    end
end
%% get
methods
    function hollow = get.hollow(BasisFunction)
        try
            hollowBFuncL = BasisFunction.BFuncL.hollow;
        catch
            hollowBFuncL = false;
        end
        if hollowBFuncL || BasisFunction.spin.hollow || isnan(BasisFunction.coe)
            hollow = true;
        else
            hollow = false;
        end
    end
    function FuncNum = get.FuncNum(BasisFunction)
        FuncNum = length(BasisFunction.BFuncL);
    end
    function BFclassL = get.BFclassL(BasisFunction)
        NumFunc = BasisFunction.FuncNum;
        BFclassL(NumFunc) = string(class(BasisFunction.BFuncL{NumFunc}));
        for i = 1:BasisFunction.FuncNum-1
            BFclassL(i) = string(class(BasisFunction.BFuncL{i}));
        end
    end
    %         function Expression = get.Expression(BasisFunction)
    %             Expression = string(BasisFunction);
    %         end
end
%% disp
methods
 dispAll(BasisFunction)
 disporbL(BasisFunction)
end
%% overload
methods
 C = eq(A,B,options)
 C = innertimes(A,B,options)
end
%% contract
methods
 B = contractrow(A,options)
end
methods(Static)
 [BFuncL] = introduce_coe(BFuncL,coeL)
 [coeL,BFuncL] = extract_coe(BFuncL,options)
end
%% rotation
methods
 U = rotation(A,Rc,Rf,tf,optionsConvection,optionsOper,optionsRm,options)
 Am = rotate(A,Rc,Rf,tf,rightorleft,optionsOper,options)
 A_Lj = rotaterow(A,Rc,Rf,tf,rightorleft,options)
 BasisFunction = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
end
methods(Static)
 [BFuncL,coeL] = rotation_func(BFuncL,R,t,options)
 spinL = rotation_spin(spin,R,t,options)
 BForb = rotation_orb(BForb,Rf,tf,options)
end
methods(Static)
end
end
