function nRemaining = vcDeleteObject(objType, val)
% Delete the specified object type of at position val
%
% Syntax:
%   nRemaining = vcDeleteObject(objType, [val])
%
% Description:
%    The ISET objects that can be deleted by this call are:
%    SCENE, OPTICALIMAGE (OI), VCIMAGE (IMGPROC), ISA (SENSOR).
%
%    If no val, this call is equivalent to vcDeleteSelectedObject
%
%    This function contains examples of usage. To access these examples,
%    type 'edit vcDeleteObject.m' into the Command Window.
%
% Inputs:
%    objType    - String. A string describing the object type.
%    val        - (Optional) Numeric. Which of the objects to remove.
%
% Outputs:
%    nRemaining - The remaining objects.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/02/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Broken Example
    vcDeleteObject('SCENE', 1);
    sceneWindow();
	vcDeleteObject('SCENE', 3);
    sceneWindow();
%}

% Get the selected object data structure and its position (val) in the list
if notDefined('objType'), error('Object type required'); end
objType = vcEquivalentObjtype(objType);

% Set val to as the selected object
if exist('val', 'var') && ~isempty(val)
    vcSetSelectedObject(objType, val)
end

% Delete the selected object
nRemaining = vcDeleteSelectedObject(objType);

end