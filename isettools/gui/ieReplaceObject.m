function ieReplaceObject(obj, val)
% Replace an object and set as selected in the vcSESSION variable
%
% Syntax:
%   ieReplaceObject(obj, [val])
%
% Description:
%    Replace an existing object in the vcSESSION global variable.
%    The object type can be SCENE, VCIMAGE, OPTICALIMAGE, or ISA.
%    The val should be the value of the object that will be replaced.
%
%    If there are no objects in the vcSESSION variable then this one becomes
%    the first entry, replacing nothing.
%
%    The function contains examples of usage below. To access these
%    examples, type 'edit ieReplaceObject.m' into the Command Window.
%
% Inputs:
%    obj - Object. The object to replace
%    val - (Optional) String. The desired object type. The default is to
%          use the current object's type.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/02/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Broken example
    ieReplaceObject(oi, 3);
    ieReplaceObject(ISA, val);
%}

if notDefined('obj'), errordlg('Object must be defined.'); end
objType = vcGetObjectType(obj);

if notDefined('val'), val = vcGetSelectedObject(objType); end
if isempty(val), val = 1; end

vcReplaceObject(obj, val);
vcSetSelectedObject(objType, val);

end