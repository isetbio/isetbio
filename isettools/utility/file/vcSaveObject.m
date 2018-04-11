function fullName = vcSaveObject(obj, fullName)
% Deprecated. Use vcExportObject. Save an ISET object into a .mat file. 
%
% Syntax:
%   fullName = vcSaveObject(obj, [fullName]);
%
% Description:
%    Users should generally use vcExportObject. If you absolutely must use
%    this routine:
%
%    The object is saved with the proper type of name so that it can be
%    loaded using vcLoadObject() at a later. When called directly, the
%    parameters and data are all saved.
%
%    If you wish to save only the parameters, without the data, then use
%    vcExportObject and set the clearDataFlag. 
%
%    The ISET objects that can be saved are: scene, opticalImage, optics, 
%        isa, pixel, or vcImage.
%
%    Examples are located within the code. To access the examples, type
%    'edit vcSaveObject.m' into the Command Window.
%
% Inputs:
%    obj      - The object you wish to save.
%    fullName - Full file name and path
%
% Outputs:
%    fullName - Full file name and path
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcExportObject
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/29/17  jnm  Formatting
%    01/22/18  dhb  Fix example so it runs.
%    01/29/18  jnm  Formatting update to match the Wiki.

% Examples:
%{
    scene = sceneCreate;
    oi = oiCreate;
    fullName = vcSaveObject(scene, fullfile(tempdir, 'myScene.mat'));
    fullName = vcSaveObject(oi, fullfile(tempdir, 'myOi.mat'))
%}

if notDefined('obj'), error('Object required.'); end

objType = vcGetObjectType(obj); 
objType = vcEquivalentObjtype(objType);

if notDefined('fullName')
    fullName = vcSelectDataFile(objType, 'w', 'mat'); 
end
if isempty(fullName), return; end

switch(lower(objType))
    case 'scene'
        scene = obj;
        save(fullName, 'scene');
    case 'optics'
        optics = obj;
        save(fullName, 'optics');
    case 'opticalimage'
        opticalimage = obj;
        save(fullName, 'opticalimage');
    case 'isa'
        sensor = obj;
        save(fullName, 'sensor');
    case 'pixel'
        pixel = obj;
        save(fullName, 'pixel');
    case 'vcimage'
        vcimage = obj;
        save(fullName, 'ip');
    otherwise
        error('Unknown object type');
end

end