function val = vcNewObjectValue(objType)
% Return an integer for a new ISET object
%
%   val = vcNewObjectValue(objType)
%
% The integer is just one more than the number of objects currently
% existing of that type.
%
% For example:
%   nextFreeValue = vcNewObjectValue('SCENE');
%
% Copyright ImagEval Consultants, LLC, 2005.

object = vcGetObjects(objType);
if isempty(object{1}),  val = 1;
else                    val = length(object) + 1;
end

return;

