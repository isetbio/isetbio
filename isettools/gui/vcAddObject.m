function val = vcAddObject(obj)
%Add and select an object to the vcSESSION data
%
%    val = vcAddObject(obj)          
%
% The object is added to the vcSESSION global variable. The object type
% can be one of the ISET object types, 
%
%   SCENE, OPTICALIMAGE, OPTICS, ISA/SENSOR, PIXEL, IP/VCI 
%
% or their aliased names in vcEquivalentObjtype
%
% The new object value is assigned the next available (new) value. 
% To see the object in the appropriate window, you call the window
% itself.
%
% Example:
%  scene = sceneCreate; 
%  newObjVal = vcAddObject(scene);
%  sceneWindow;
%
% See also:  vcAddAndSelectObject.m (
%  
% Copyright ImagEval Consultants, LLC, 2013

%%
global vcSESSION;

% Get a value
% Makes objType proper type and forces upper case.
objType = obj.type;
val = vcNewObjectValue(objType);

% gather to avoid distributed component (e.g. gpuArray)
obj = gatherStruct(obj);

%% Assign object to the vcSESSION global.

% Should be ieSessionSet, not this.
if exist('obj','var')
    switch lower(objType)
        case {'scene'}
            lum = sceneGet(obj, 'luminance');
            obj = sceneSet(obj, 'luminance', lum);
            vcSESSION.SCENE{val} = obj;
        case {'opticalimage'}
            vcSESSION.OPTICALIMAGE{val} = obj;
        case {'optics'}
            oi = vcSESSION.OPTICALIMAGE{val};
            oi = oiSet(oi,'optics',obj);
            vcSESSION.OPTICALIMAGE{val} = oi;
        case {'sensor'}
            vcSESSION.ISA{val} = obj;
        case {'pixel'}
            sensor = vcSESSION.ISA{val};
            sensor = sensorSet(sensor,'pixel',obj);
            vcSESSION.ISA{val} = sensor;
        case {'vcimage'}
            vcSESSION.VCIMAGE{val} = obj;
        case {'display'}
            vcSESSION.DISPLAY{val} = obj;
        otherwise
            error('Unknown object type');
    end
end

vcSetSelectedObject(objType,val);


return;
