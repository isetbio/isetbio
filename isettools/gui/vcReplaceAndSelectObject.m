function vcReplaceAndSelectObject(obj, val)
% Replace an object and set as selected in the vcSESSION variable
%
% Syntax:
%   vcReplaceAndSelectObject(obj, [val])
%
% Description:
%    Replace an existing object in the vcSESSION global variable.
%    The object type can be SCENE, VCIMAGE, OPTICALIMAGE, or ISA.
%    The val should be the value of the object that will be replaced.
%
%    If there are no objects in the vcSESSION variable then this one
%    becomes the first entry, replacing nothing.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcReplaceAndSelectObject.m' into the Command Window.
%
% Inputs:
%    obj - Object. The desired object.
%    val - (Optional) Numeric. The index of the object to be replaced.
%          Default is 1.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/11/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - skipping broken example
    vcReplaceAndSelectObject(oi, 3);
    vcReplaceAndSelectObject(ISA, val);
%}

if notDefined('obj'), errordlg('Object must be defined.'); end
objType = vcGetObjectType(obj);

if notDefined('val'), val = vcGetSelectedObject(objType); end
if isempty(val), val = 1; end

vcReplaceObject(obj, val);
vcSetSelectedObject(objType, val);

end