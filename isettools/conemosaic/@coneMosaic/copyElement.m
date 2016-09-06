
function cpObj = copyElement(obj)
% make a shallow copy of the properties
cpObj = copyElement@matlab.mixin.Copyable(obj);

% make deep copy of the cone and macular class
cpObj.pigment = cpObj.pigment.copy();
cpObj.macular = cpObj.macular.copy();
end