function [val, sOBJECT] = vcGetSelectedObject(objType)
% Returns the number of the currently selected object of some type
%
% Syntax:
%   [val, sOBJECT] = vcGetSelectedObject(objType)
%
% Description:
%    This routine returns both the number and object itself, if requested.
%    You may wish to use vcGetObject(objType) in general. That routine
%    returns the currently selected object, of objType.
%
%    The set of possible objects returned are: scene, pixel, optics,
%    {opticalImage, oi}, vcImage, graphWin, {isa, sensor}
%
%    Originally, this routine was the main one used to return objects.
%
%        val = vcGetSelectedObject('SCENE')
%        [val, pixel] = vcGetSelectedObject('PIXEL')
%        [val, sensor] = vcGetSelectedObject('SENSOR')
%        [val, vci] = vcGetSelectedObject('VCIMAGE')
%        [val, vci] = vcGetSelectedObject('IMGPROC')
%        [val, display] = vcGetSelectedObject('DISPLAY')
%
%    As of May 2004, I started using
%       obj = vcGetObject(objType)
%       obj = vcGetObject(objType, val)
%
%    Currently, I only call this routine directly if I need the val
%
%       val = vcGetSelectedObject('foo')
%
% Inputs:
%    objType - String. A string describing the object type.
%
% Outputs:
%    val     - Numeric. The index of the object of the specified type.
%    sObject - Object. The desired object.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcGetObject

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/11/18  jnm  Formatting

%%
global vcSESSION

objType = vcEquivalentObjtype(objType);
val = [];

switch lower(objType)
    case 'scene'
        if checkfields(vcSESSION, 'SELECTED', 'SCENE')
            val = vcSESSION.SELECTED.SCENE;
        end
    case {'opticalimage', 'oi'}
        if checkfields(vcSESSION, 'SELECTED', 'OPTICALIMAGE')
            val = vcSESSION.SELECTED.OPTICALIMAGE;
        end
    case {'optics'}
        if checkfields(vcSESSION, 'SELECTED', 'OPTICALIMAGE')
            val = vcSESSION.SELECTED.OPTICALIMAGE;
        end
    case {'isa', 'sensor'}
        if checkfields(vcSESSION, 'SELECTED', 'ISA')
            val = vcSESSION.SELECTED.ISA;
        end
    case {'vcimage', 'ip', 'imgproc'}
        if checkfields(vcSESSION, 'SELECTED', 'VCIMAGE')
            val = vcSESSION.SELECTED.VCIMAGE;
        end
    case {'pixel'}
        if checkfields(vcSESSION, 'SELECTED', 'ISA')
            val = vcSESSION.SELECTED.ISA;
        end
    case {'display'}
        if checkfields(vcSESSION, 'SELECTED', 'DISPLAY')
            val = vcSESSION.SELECTED.DISPLAY;
        end
    otherwise
        error('Unknown object type.');
end

if nargout == 2, sOBJECT = vcGetObject(objType, val); end

end