function [sOBJECT, val] = vcGetObject(objType, val)
% Retrieve an object from vcSESSION structure
%
% Syntax:
%   [sOBJECT, val] = vcGetObject(objType, [val])
%
% Description:
%   Return an ISETBIO object from the database.
%
%    This routine replaces: [val, sOBJECT] = vcGetSelectedObject('SCENE');
%
%    If you need the val of the currently selected object, you can use
%
%    [obj, val] = vcGetObject('SCENE');
%
% Inputs:
%    objType - {SCENE, PIXEL, OPTICS, OI, SENSOR, IP, GRAPHWIN}
%    val     - Which object in the list.  If not passed in, the
%              currently 'selected' object is returned.
%
% Outputs:
%    sObject - String. String describing object type
%    val     - Object. The object in question.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%  vcGetObjects, vcGetROIData
%

%%
global vcSESSION

% For speed
if ~exist('objType', 'var') || isempty(objType)
    error('objType must be defined');
end
if ~exist('val', 'var') || isempty(val)
    val = vcGetSelectedObject(objType);
end

% Translates the common names (e.g., sensor or ip) into the official name.
objType = vcEquivalentObjtype(objType);

%%
if ~isempty(val)
    switch(lower(objType))
        case {'scene', 'isa', 'opticalimage', 'vcimage', 'display'}
            eval(['sOBJECT = vcSESSION.', objType, '{val};']);
        case {'optics'}
            sOBJECT = oiGet(vcSESSION.OPTICALIMAGE{val}, 'optics');
        otherwise
            error('Unknown object type.');
    end
else
    % No val.  Return empty.
    sOBJECT = [];
end

end
