function val = vcNewObjectValue(objType)
% Return an integer for a new ISET object
%
% Syntax:
%   val = vcNewObjectValue(objType)
%
% Description:
%    The integer is just one more than the number of objects currently
%    existing of that type.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcNewObjectValue.m' into the Command Window.
%
% Inputs:
%    objType - String. A string describing the object's type.
%
% Outputs:
%    val     - Integer. An integer describing the number of objects of the
%              specified type.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/09/18  jnm  Formatting

% Examples:
%{
    nextFreeValue = vcNewObjectValue('SCENE');
%}

object = vcGetObjects(objType);

% It seems that ieInit was not run.  So return and tell the user.
if isempty(object)
    val = [];
    return;
end

if isempty(object{1}), val = 1; else, val = length(object) + 1; end

return;
