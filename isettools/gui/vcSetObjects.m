function vcSetObjects(objType, newObj)
% Set the cell array of a particular type of object.
%
% Syntax:
%   vcSetObjects(objType, newObj);
%
% Description:
%    Set the cell array of a particular type of object.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcSetObjects.m' into the Command Window.
%
% Inputs:
%    objType - String. A string describing the object's type.
%    newObj  - Object. The object to replace with.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Examples:
%{
    % ETTBSkip - Skipping broken example
    scene = vcSESSION.SCENE;
    n = length(scene);
    scene = scene{1:n - 1};
    vcSetObject('SCENE', scene);
%}

global vcSESSION;

objType = vcEquivalentObjtype(objType);

eval(['vcSESSION.', objType, ' = newObj;']);

return;
