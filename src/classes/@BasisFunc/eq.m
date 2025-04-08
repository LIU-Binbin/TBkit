function C = eq(A,B,options)
arguments
A BasisFunc;
B BasisFunc;
options.strict = false;
options.spin = false;
options.orb = true;
options.BFuncL =true;
end
C = true;
if class(A.BFuncL) ~= class(B.BFuncL)
C =false;
return;
end
if options.spin
if ~eq(A.spin,B.spin,'strict',options.strict)
C =false;
return;
end
end
if options.orb
if isempty(A.BForb) &&  ~isempty(B.BForb)
C =false;
return;
end
if ~isempty(A.BForb) &&  isempty(B.BForb)
C =false;
return;
end
if ~Oper.allclose(A.BForb,B.BForb)
C =false;
return;
end
end
if options.BFuncL
if ~eq(A.BFuncL,B.BFuncL,'strict',options.strict)
C =false;
return;
end
end
end
