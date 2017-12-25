function [newVal, fullName] = vcImportObject(objType, fullName, ...
    preserveDataFlag)
% Import an ISET structure from a file to the vcSESSION structure 
%
% Syntax:
%   [newVal, fullFileName] = vcImportObject([objType], [fullName], ...
%        [preserveDataFlag])
%
% Description:
%    The parameters of an ISET object, saved in a file, are read and
%    attached to the variable, vcSession. 
%
%    Imported objects do not have any image data attached to them. If a
%    pixel or optics are imported, the data in the current ISA or OI are
%    cleared and must be recomputed. This is done to assure consistency
%    between the data and the structure parameters. 
%
%    This function also loads optics, attaching them to the current optical
%    image. The default for optics is to preserve the data. The default for
%    other objects is to clear the data.
%
% Inputs:
%    objType          - (Optional) Object type to import. Default is scene.
%    fullName         - (Optional) Full path and file name. Default is []
%    preserveDataFlag - (Optional) Whether or not to preserve prior data.
%                       Default depends on the object type. 1 for optics
%                       and pixel cases, and 0 for everything else.
%
% Outputs:
%    newVal           - The new value for the variable
%    fullName         - the full file name and path.
%
% Notes:
%    * [Note: JNM - TODO: unify vcLoad (below) and vcLoadObject, as per
%      below note in the function.]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/28/17  jnm  Formatting
%    11/29/17  jnm  Note about example

% Examples:
%{
    scene = sceneCreate;
    fname = tempname;
    save(fname,'scene');
    newVal = vcImportObject('SCENE',fname);
%}

if notDefined('objType'), objType = 'SCENE'; end
if notDefined('fullName'), fullName = []; end
if notDefined('preserveDataFlag')
    switch lower(objType)
        case {'scene', 'opticalimage', 'oi', 'isa', 'vcimage'}
            preserveDataFlag = 0;
        otherwise
            % optics and pixel case, but I don't think the pixel has data.
            preserveDataFlag = 1;
    end
end

% [Note: XXX - that there is vcLoad in this file and vcLoadObject is a
% different function. Should be unified some day, sigh.]
switch lower(objType)
    case {'scene', 'opticalimage', 'oi', 'isa', 'vcimage'}
        % Load the object into a new value assigned by vcLoadObject.
        [newVal, fullName] = vcLoadObject(objType, fullName);
        if isempty(newVal), return; end
    case {'optics'}
        [newVal, oi] = vcGetSelectedObject('OPTICALIMAGE');
        [optics, fullName] = vcLoad(objType, fullName);
        if ~isempty(optics)
            oi = oiSet(oi, 'optics', optics);
            if ~preserveDataFlag
                oi = oiClearData(oi);
            end
            vcReplaceAndSelectObject(oi, newVal);
        end
    otherwise
        error('Unknown object type.');
end


end

%----------------------------------------------------------
function [obj, fullName] = vcLoad(objType, fullName)
% Function to handle loading pixels, optics, and (in the future) displays.
%
% Syntax:
%   [obj, fullName] = vcLoad(objType, fullName)
%
% Description:
%    This routine will handle loading pixels, optics, and in the future, it
%    will handle displays.
%
% Inputs:
%    objType  - The required object type
%    fullName - the full file name and path
%
% Outputs:
%    obj      - The object in question
%    fullName - The object path and filename
%

obj = [];

if notDefined('fullName')
    windowTitle = sprintf('Select %s file name', objType);
    fullName = vcSelectDataFile('session', 'r', 'mat', windowTitle);
    if isempty(fullName), return; end
end

switch(lower(objType))      
    case 'optics'
        data = load(fullName, 'optics');
        obj = data.optics;
        
    otherwise
        error('Unknown object type');
end

end