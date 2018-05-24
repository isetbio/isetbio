function val = vcAddObject(obj)
% Add and select an object to the vcSESSION data
%
% Syntax:
%   val = vcAddObject(obj)
%
% Description:
%    The object is added to the vcSESSION global variable. The object type
%    can be one of the ISET object types,
%
%      SCENE, OPTICALIMAGE, OPTICS, ISA/SENSOR, PIXEL, IP/VCI
%
%    or their aliased names in vcEquivalentObjtype
%
%    The new object value is assigned the next available (new) value. To
%    see the object in the appropriate window, you call the window itself.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcAddObject.m' into the Command Window.
%
% Inputs:
%    obj - Object. The object in question.
%
% Outputs:
%    val - Numeric. The index of the object.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcAddAndSelectObject.m
%

% History:
%    xx/xx/13       Copyright ImagEval Consultants, LLC, 2013
%    05/10/18  jnm  Formatting

% Examples:
%{
    scene = sceneCreate;
    newObjVal = vcAddObject(scene);
    sceneWindow;
%}

%%
global vcSESSION;

% Get a value
% Makes objType proper type and forces lowercase.
objType = obj.type;
val = vcNewObjectValue(objType);

% gather to avoid distributed component (e.g. gpuArray)
obj = gatherStruct(obj);

%% Assign object to the vcSESSION global.
% Should be ieSessionSet, not this.
if exist('obj', 'var')
    switch lower(objType)
        case {'scene'}
            lum = sceneGet(obj, 'luminance');
            obj = sceneSet(obj, 'luminance', lum);
            vcSESSION.SCENE{val} = obj;
        case {'opticalimage'}
            vcSESSION.OPTICALIMAGE{val} = obj;
        case {'optics'}
            oi = vcSESSION.OPTICALIMAGE{val};
            oi = oiSet(oi, 'optics', obj);
            vcSESSION.OPTICALIMAGE{val} = oi;
        case {'sensor'}
            vcSESSION.ISA{val} = obj;
        case {'vcimage'}
            vcSESSION.VCIMAGE{val} = obj;
        case {'display'}
            vcSESSION.DISPLAY{val} = obj;
        otherwise
            error('Unknown object type');
    end
end

vcSetSelectedObject(objType, val);

return;
