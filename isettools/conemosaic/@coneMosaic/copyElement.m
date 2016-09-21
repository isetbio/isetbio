function cpObj = copyElement(obj)
% Make a shallow copy of the properties
% 
% HJ ISETBIO Team 2016

    cpObj = copyElement@matlab.mixin.Copyable(obj);
    
    % make deep copy of the cone and macular class
    cpObj.pigment = cpObj.pigment.copy();
    cpObj.macular = cpObj.macular.copy();
end