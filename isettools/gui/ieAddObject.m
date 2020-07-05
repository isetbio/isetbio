function val = ieAddObject(obj)
% Add and select an object to the vcSESSION data
%
% Syntax:
%   val = ieAddObject(obj)
%
% Description:
%    The object is added to the vcSESSION global variable. The object type
%    can be one of the ISET object types,
%       - SCENE, OPTICALIMAGE, OPTICS, ISA/SENSOR, PIXEL, IP/VCI
%       - Or their aliased names in vcEquivalentObjtype
%
%    The new object value is assigned the next available (new) value. To
%    see the object in the appropriate window, you call the window itself.
%
%    There are examples contained in the code. To access the examples, type
%    'edit ieAddObject.m' into the Command Window.
%
% Inputs:
%    obj - Object. The object to add to the session data
%
% Outputs:
%    val - The value within the session that the aforementioned object has
%          been assigned to.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Determine if "should be ieSessionSet, not this." changes need
%      to be implemented.
%
% See Also:
%    vcAddAndSelectObject.m
%

% History:
%    xx/xx/13       Copyright ImagEval Consultants, LLC, 2013
%    02/28/18  jnm  Formatting

% Examples:
%{
    scene = sceneCreate;
    newObjVal = ieAddObject(scene);
    sceneWindow;
%}
%%
global vcSESSION;

% Get a value
% Makes objType proper type and forces upper case.
objType = obj.type;
val = vcNewObjectValue(objType);
if isempty(val), disp('Please initialize ISET (ieInit)'); return; end

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
        case {'vcimage'}
            vcSESSION.VCIMAGE{val} = obj;
        case {'display'}
            vcSESSION.DISPLAY{val} = obj;
        otherwise
            error('Unknown object type');
    end
end

vcSetSelectedObject(objType, val);

end
