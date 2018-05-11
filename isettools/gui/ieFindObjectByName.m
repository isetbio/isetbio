function val = ieFindObjectByName(objType, objName)
% Return an object number given its type and name
%
% Syntax:
%   val = ieFindObjectByName(objType, objName)
%
% Description:
%    Find the object of objType with name of objName. Return its slot in
%    the cell array within vcSESSION.
%
%    The code below contains examples of function usage. To access, type
%    'edit ieFindObjectByName.m' into the Command Window.
%
% Inputs:
%    objType - The type of the object
%    objName - The name of the object
%
% Outputs:
%    val     - The number of the object
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Function seems .... tempermental? See second example.]

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/28/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - no vcimage in existence prior to this call
    val = ieFindObjectByName('vcimage', 'mono1')
%}
%{
    s1 = sceneCreate();
    s2 = sceneCreate('gridlines');
    val = ieFindObjectByName('scene', 'gridlines')
    % [Note: JNM - val should be '2', but is []???]
%}

obj = vcGetObjects(objType);
if length(obj) == 1 && isempty(obj{1})
    warning('No objects of type %s', objType);
    return;
end

val = [];
for ii = 1:length(obj)
    if strcmp(obj{ii}.name, objName)
        val = ii;
        return;
    end
end

return;