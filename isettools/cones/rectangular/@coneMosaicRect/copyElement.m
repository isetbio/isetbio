function cpObj = copyElement(obj)
% Make a shallow copy of the properties
%
% Syntax:
%  cpObj = copyElement(obj)
%
% Description:
%    Copy an element out of an object?
%
% Inputs:
%    obj   - The cone mosaic object you wish to copy
%
% Outputs:
%    cpObj - The created copy of the cone mosaic object
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: DHB - WHAT IS A SHALLOW COPY? WHY DO THE COMMENTS IN THE
%     ACTUAL CODE SAY IT IS MAKING A DEEP COPY? WHEN WOULD ONE USE THIS
%     ROUTINE, AND WHY?]
%    * [Note: JNM - I think that you create a shallow copy of the general
%      object and a deep copy of the particular classes listed below? But I
%      am definitely unsure about that.]

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/23/18  jnm  Formatting

    cpObj = copyElement@matlab.mixin.Copyable(obj);

    % Make deep copy of the cone and macular class
    cpObj.pigment = cpObj.pigment.copy();
    cpObj.macular = cpObj.macular.copy();
end