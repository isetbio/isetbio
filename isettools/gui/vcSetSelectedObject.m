function vcSetSelectedObject(objType, val)
% Set selected object in vcSESSION variable.
%
% Syntax:
%   vcSetSelectedObject(objType, val)
%
% Description:
%    Sets the currently selected object, where the ISET object type might
%    be SCENE, OPTICALIMAGE, ISA, and VCIMAGE.
%
%    If val is 0 (or less than one) the selected value is set to empty.  If
%    it is a positive integer, then we check that it is in range of the
%    number of objects of that type, warn if it isn't, and go ahead and set
%    the value.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcSetSelectedObject.m' into the Command Window.
%
% Inputs:
%    objType - String. A string describing the object type.
%    val     - Numeric. A value to indicate the object in the list of its
%              type. If < 1, then this is empty. If the number doesn't
%              exist, create the object at that value.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: We should eliminate these routines and instead do the work
%      through switch statements that go to ieSessionSet
%
% See Also:
%    vcAddAndSelectObject, vcNewObjectValue
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/10/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Skip broken example
    vcSetSelectedObject('SCENE', 1)
    vcSetSelectedObject('OPTICALIMAGE', 3);
    vcSetSelectedObject('OPTICALIMAGE', 3);
%}

global vcSESSION;

objType = vcEquivalentObjtype(objType);

if isempty(val) || val < 1
    eval(['vcSESSION.SELECTED.', objType, '= [];'])
else
    nObjects = length(vcGetObjects(objType));
    if val <= nObjects
        eval(['vcSESSION.SELECTED.', objType, '= val;']);
    else
        error('Selected object out of range.');
    end
end

return;
