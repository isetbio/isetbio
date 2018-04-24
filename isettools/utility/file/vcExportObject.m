function fullName = vcExportObject(obj, fullName, clearDataFlag)
% Save an object into the relevant object directory
%
% Syntax:
%   fullName = vcExportObject(obj, [fullName], [clearDataFlag=0]);
%
% Description:
%    Save the parameters of a vcSession object in a .mat file. The object
%    parameters can be loaded at a later time. In the future, we will
%    support exporting to formats other than Matlab files, which is why
%    this routine exists. And we may allow exporting the data as well.
%
%    If fullName is not specified, a GUI opens to choose the name.
%
%    Examples are located within the code. To access the examples, type
%    'edit vcExportObject.m' into the Command Window.
%
% Inputs:
%    obj           - object to add to the file
%    fullName      - (Optional) full file path and name. Default is blank.
%    clearDataFlag - (Optional) Boolean whether or not to clear the data
%                    after exporting it. Default is 0 (false).
%
% Outputs:
%    fullName      - full file path and name of the file where the object
%                    has been saved.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/28/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match the Wiki.
  
% Examples:
%{
    scene = sceneCreate;
    fullName = vcExportObject(scene, fullfile(tempdir, 'myCompany'));
%}

if notDefined('obj'), error('You must define an object to save.'); end
if notDefined('clearDataFlag'), clearDataFlag = 0;   end
if notDefined('fullName'), fullName = []; end

objType = obj.type;
objType = vcEquivalentObjtype(objType);

switch(lower(objType))     
    case {'scene', 'opticalimage', 'isa', 'vcimage'}
        if clearDataFlag, obj.data = []; end
    case 'optics'
        error('NYI');
    otherwise
        error('Unknown object type');
end

fullName = vcSaveObject(obj, fullName);

end