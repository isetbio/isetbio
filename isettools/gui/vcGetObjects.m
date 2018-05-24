function [obj, nObj] = vcGetObjects(objType)
% Retrieve cell larray of objects of objType
%
% Syntax:
%   [obj, nObj] = vcGetObjects(objType)
%
% Description:
%    Return the cell array of a particular type of object.  Optionally, the
%    total number of objects of this type is returned, too.
%
%    The code contains examples of usage. To access, type 'edit
%    vcGetObjects.m' into the Command Window.
%
% Inputs:
%    objType - String. A string describing the object type.
%
% Outputs:
%    obj     - Cell Array. Cell array of objects.
%    nObj    - Numeric. The number of objects in the cell array.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcGetObject, vcGetSelectedObject
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/03/18  jnm  Formatting

% Examples:
%{
    obj = vcGetObjects('PIXEL');
    [obj, nObj] = vcGetObjects('SCENE');
%}

global vcSESSION;

% Translate various names into the proper name used in vcSESSION.
% This routine also forces upper case on the object type, as required.
objType = vcEquivalentObjtype(objType);

if checkfields(vcSESSION, objType)
    eval(['obj = vcSESSION.', objType, ';']);
else
    obj{1} = [];
end

if nargout == 2, nObj = length(obj); end

return