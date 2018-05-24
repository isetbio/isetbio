function nObj = vcCountObjects(objType)
% Count how many objects of a type are contained in vcSESSION
%
% Syntax:
%   nObj = vcCountObjects(objType)
%
% Description:
%    This routine returns the number of objects of objType in the vcSESSION
%    structure. If the structure is empty, 0 is return.
%
% Inputs:
%    objType - String. The object type.
%
% Outputs:
%    nObj    - Numeric. The number of objects of that type that exist.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/02/18  jnm  Formatting.

obj = vcGetObjects(objType);
if isempty(obj{1}), nObj = 0; else, nObj = length(obj); end
return;
