function [newVal, fullName] = vcLoadObject(objType, fullName, val)
% Deprecated. Use vcImportObject. Load ISET object into vcSESSION structure
%
% Syntax:
%   [newVal, fullName] = vcLoadObject(objType, [fullName], [val])
%
% Description:
%    Users should use vcImportObject instead. That routine calls this one
%    when appropriate, but handles more types of objects.
%
%    If you must use this routine:
%
%    ISET objects can be saved as Matlab (.mat) files and then loaded. This
%    routine imports the data from an ISET object.
%
%    The data are loaded and stored in the vcSession data structure. You
%    can assign it a particular slot in the cell array of objects using
%    val, or if no val is passed to this routine the object will be given a
%    new value. The file name used to import the data is also assigned to
%    the object name.
%
%    The object types that can be loaded are scene, opticalImage, isa, or
%    vcImage
%
% Inputs:
%    objType  - (Optional) The object type. Default is 'session'.
%    fullName - (Optional) Full file name and path. Default is to query the
%               user with a GUI.
%    val      - (Optional) Slot in the cell array of objects in which to
%               save the data. Default is a new value assigned in-function. 
%
% Outputs:
%    newVal   - The val, if provided, else a new value, indicating which
%               slot in the vcSession data structure.
%    fullName - Full file name and path
%
% Notes:
%    * [Note: JNM - Should we perhaps use a more concrete example, rather
%      than querying the user to select a file?]
%    * [Note: JNM - Was the 'isa' problem resolved? Is the variable naming
%      problem related to a function with the same name?]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/30/17  jnm  Formatting

% Examples:
%{
    newVal = vcLoadObject('scene');
%}

%%
if notDefined('objType'), objType = 'scene'; end
if notDefined('val')
    newVal = vcNewObjectValue(objType);
else
    newVal = val;
end

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
        % We will have to make many changes to fix this. In  Matlab 7.04
        % the load of isa generates a warning, and the Mathworks kindly
        % changes my variable name. We trap this condition here and handle
        % it. But if we have to do it lots of places, we are in trouble.
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