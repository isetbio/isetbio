function cpObj = copyElement(obj)
%COPYELEMENT  Make a shallow copy of the properties
%
% [DHB NOTE: WHAT IS A SHALLOW COPY? WHY DO THE COMMENTS IN THE ACTUAL CODE SAY
% IT IS MAKING A DEEP COPY? WHEN WOULD ONE USE THIS ROUTINE, AND WHY?]

% HJ ISETBIO Team 2016

    cpObj = copyElement@matlab.mixin.Copyable(obj);
    
    % Make deep copy of the cone and macular class
    cpObj.pigment = cpObj.pigment.copy();
    cpObj.macular = cpObj.macular.copy();
end