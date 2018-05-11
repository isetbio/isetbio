function objType = vcGetObjectType(obj)
% Return a string describing the type of ISET object.
%
% Syntax:
%   objType = vcGetObjectType(obj)
%
% Description:
%    The set of types includes:  SCENE, OPTICALIMAGE, ISA, VCIMAGE, PIXEL
%    and OPTICS. This routine will perform the task properly for any
%    structure that as a .type field.
%
%    The code contains examples of function usage. To access these
%    examples, type 'edit vcGetObjectType.m' into the Command Window.
%
% Inputs:
%    obj     - Object. The object in question.
%
% Outputs:
%    objType - String. A string describing the object type.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/03/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - skip broken examples.
    oi = vcGetObject('oi');
    vcGetObjectType(oi)

    obj = oiGet(oi, 'optics');
    vcGetObjectType(obj)
%}

if checkfields(obj, 'type')
    objType = obj.type;
    return;
else
    error('Object does not have a type field.');
end
