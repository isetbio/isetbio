function val = vcAddAndSelectObject(objType, obj)
% Add and select an object to the vcSESSION data
%
% Syntax:
%   val = vcAddAndSelectObject(obj)           % (preferred)
%   val = vcAddAndSelectObject(objType, obj)  % (supported)
%
% Description:
%    The object is added to the vcSESSION global variable. The object type
%    can be one of the ISET object types, SCENE, VCIMAGE, OPTICALIMAGE
%    ISA/SENSOR or their equivalents.
%
%    The new object value is assigned the next new value. To see the object
%    in the appropriate window, you can call the window itself.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcAddAndSelectObject.m' into the Command Window.
%
% Inputs:
%    objType - String. A string describing the object's type.
%    obj     - Object. The object in question.
%
% Outputs:
%    val     - Numeric. The instance of the object of its type.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Start using ieSessionSet instead of direct assignments here.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/10/18  jnm  Formatting

% Examples:
%{
    scene = sceneCreate;
    newObjVal = vcAddAndSelectObject(scene);
    sceneWindow;
%}
%{
    % Older syntax is supported, but not preferred
    oi = oiCreate;
    scene = sceneCreate;
    newObjVal = vcAddAndSelectObject('OPTICALIMAGE', oi);
    vcAddAndSelectObject('SCENE', scene);
%}

global vcSESSION;

% This way, the call can be
% vcAddAndSelectObject(scene) instead of
% vcAddAndSelectObject('scene', scene), which was the original
%
if checkfields(objType, 'type'), obj = objType; objType = objType.type; end

% gather to avoid distributed component (e.g. gpuArray)
obj = gatherStruct(obj);

% Makes objType proper type and forces upper case.
objType = vcEquivalentObjtype(objType);
val = vcNewObjectValue(objType);

% Assign object, passed in as 3rd variable, to the vcSESSION global.
% Should be ieSessionSet, not this.
if exist('obj', 'var')
    switch upper(objType)
        case {'SCENE'}
            vcSESSION.SCENE{val} = obj;
        case {'OPTICALIMAGE'}
            vcSESSION.OPTICALIMAGE{val} = obj;
        case {'ISA'}
            vcSESSION.ISA{val} = obj;
        case {'VCIMAGE'}
            vcSESSION.VCIMAGE{val} = obj;
        case {'DISPLAY'}
            vcSESSION.DISPLAY{val} = obj;
        otherwise
            error('Unknown object type');
    end
end

vcSetSelectedObject(objType, val);

return;
