function vcReplaceObject(obj, val)
% Replace an object in the vcSESSION variable
%
% Syntax:
%   vcReplaceObject(obj, [val])
%
% Description:
%    Replace an existing object, either a SCENE, VCIMAGE, OPTICS, PIXEL, 
%    OPTICALIMAGE, or ISA, in the vcSESSION global variable.
%
%    When  replacing OPTICS or PIXEL the val refers to the OPTICALIMAGE or
%    SENSOR that contain the OPTICS or PIXEL.
%
%    The object that is replaced (or its parent) are then set to be the
%    selected object.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcReplaceObject.m' into the Command Window.
%
% Inputs:
%    obj - Object. The object in question.
%    val - Numeric. The number of the object to be replaced. If val is not
%          specified, then the currently selected object is replaced.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/11/18  jnm  Formatting.

% Examples:
%{
    % ETTBSkip - Skipping broken examples.
    vcReplaceObject(oi, 3);
    vcReplaceObject(oi);
    vcReplaceObject(ISA, val);
%}

%%
global vcSESSION;

objType = vcGetObjectType(obj);
objType = vcEquivalentObjtype(objType);

%%
if notDefined('val')
    val = vcGetSelectedObject(objType);
    if isempty(val), val = 1; end
end

% Should be handled by ieSessionSet
switch lower(objType)
    case 'scene'
        vcSESSION.SCENE{val} = obj;
    case 'opticalimage'
        vcSESSION.OPTICALIMAGE{val} = obj;
    case 'optics'
        vcSESSION.OPTICALIMAGE{val}.optics = obj;
    case 'isa'
        vcSESSION.ISA{val} = obj;
    case 'pixel'
        vcSESSION.ISA{val}.pixel = obj;
    case 'vcimage'
        vcSESSION.VCIMAGE{val} = obj;
    otherwise
        error('Unknown object type');
end

vcSetSelectedObject(objType, val)

return;
