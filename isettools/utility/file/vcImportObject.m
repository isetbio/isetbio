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
%    Examples are located within the code. To access the examples, type
%    'edit vcImportObject.m' into the Command Window.
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
% Optional key/value pairs:
%    None.
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
    save(fname, 'scene');
    newVal = vcImportObject('SCENE', fname);
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
        [optics, fullName] = vcLoadOptics(objType, fullName);
        if ~isempty(optics)
            oi = oiSet(oi, 'optics', optics);
            if ~preserveDataFlag, oi = oiClearData(oi); end
            ieReplaceObject(oi, newVal);
        end
    otherwise
        error('Unknown object type.');
end

end

%----------------------------------------------------------
function [obj, fullName] = vcLoadOptics(objType, fullName)
% Function to handle loading pixels, optics, and (in the future) displays.
%
% Syntax:
%   [obj, fullName] = vcLoadOptics(objType, [fullName])
%
% Description:
%    This routine will handle loading pixels, optics, and in the future, it
%    will handle displays.
%
% Inputs:
%    objType  - The required object type
%    fullName - (Optional) the full file name and path. Default is to
%               query the user to select a file.
%
% Outputs:
%    obj      - The object in question
%    fullName - The object path and filename
%
% Optional key/value pairs:
%    None.
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

% To be deprecated or at least heavily changed. Loads an ISET object
% from a file into the vcSESSION structure. Rather like a read and
% ieAddObject() call.
function [newVal, fullName] = vcLoadObject(objType, fullName, val)
% Use vcImportObject, which will call this one as necessary.
%
% Syntax:
%   [newVal, fullName] = vcLoadObject([objType], [fullName], [val])
%
% Description:
%    Use vcImportObject instead. vcImportObject will call this one when it
%    is appropriate, but will handle more types of objects and a more
%    general case.
%
%    ISETBIO objects can be saved as Matlab (.mat) files and then loaded.
%    This routine imports the data from an ISET object.
%
%    The data are loaded and stored in the vcSESSION data structure. You
%    can assign it a particular slot in the cell array of objects using
%    VAL, or if no VAL is passed to this routine the object will be given a
%    new value. The file name used to import the data is also assigned to
%    the object name.
%
%    The object types that can be loaded are SCENE, OPTICALIMAGE, ISA, or
%    VCIMAGE
%
%    Examples are located within the code. To access the examples, type
%    'edit vcLoadObject.m' into the Command Window.
%
% Inputs:
%    objType  - (Optional) The object type. Default is scene.
%    fullname - (Optional) The full file name. Default is to query user.
%    val      - (Optional) The object value(s). Default is to create a new
%               value for the object.
%
% Outputs:
%    newVal   - The newly created value.
%    fullName - The full file name.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/29/18  jnm  Formatting

% Examples:
%{
    % [Note: DHB - This does not seem to actually work in that selecting
    % any of the available files in the dialog produces a warning that
    % "Variable 'scene' not found."]
    scene = sceneCreate;
    ieAddObject(scene);
    vc
    newVal = vcLoadObject('SCENE');
%}

if notDefined('objType'), objType = 'scene'; end
if notDefined('val'), newVal = vcNewObjectValue(objType);
else, newVal = val; end

% Parse the object type
objType  = vcEquivalentObjtype(objType);

% Set up the full file name
if notDefined('fullName')
    fullName = vcSelectDataFile('stayput', 'r', 'mat');
    if isempty(fullName), newVal = []; return; end
end
[~, objName] = fileparts(fullName);

%%
switch(lower(objType))
    case 'scene'
        data = load(fullName, 'scene');
        data.scene.name = objName;
        vcAddAndSelectObject('scene', data.scene);
        
    case 'opticalimage'
        data = load(fullName, 'opticalimage');
        data.opticalimage.name = objName;
        vcAddAndSelectObject('opticalimage', data.opticalimage);
        
    case 'isa'
        % We have problems with the variable name isa in 7.04. 
        % We will have to make many changes to fix this.
        % In  Matlab 7.04 the load of isa generates a warning, and the
        % Mathworks kindly changes my variable name. We trap this condition
        % here and handle it. But if we have to do it lots of places, we
        % are in trouble.
        % warning('off');
        data = load(fullName, 'isa');
        % warning('on');

        % This is what they rename the variable in Matlab 7.04
        if checkfields(data, 'isa_'), data.isa = data.isa_; end
        data.isa.name = objName;
        vcAddAndSelectObject('isa', data.isa);
        
    case 'vcimage'
        data = load(fullName, 'vcimage');
        data.vcimage.name = objName;
        vcAddAndSelectObject('vcimage', data.display);
        
    otherwise
        error('Unknown object type');
end

end