function TrueOrFalse = LatticeVectorTest(Vector,Accuracy)
if nargin < 2
    Accuracy = 1e-6;
end
TrueOrFalse = true;
for i = 1:numel(Vector)
    if abs(rem(Vector(i),1)) > Accuracy
        TrueOrFalse = false;
        return;
    end
end
end