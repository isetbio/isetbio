function objNames = vcGetObjectNames(objType)
% Compile a list of object names from vcSESSION variable
%
% Syntax:
%   objNames = vcGetObjectNames(objType)
%
% Description:
%    Compile a a list of object names from the vcSESSION variable.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcGetObjectNames.m' into the Command Window.
%
% Inputs:
%    objType  - String. A string describing the object type.
%
% Outputs:
%    objNames - List. The list of object names.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/03/18  jnm  Formatting

% Examples:
%{
    vcGetObjectNames('oi')
%}

if notDefined('objType'), objType = 'scene'; end

objects = vcGetObjects(objType);
nObj = length(objects);

% There may be just one, empty object. In which case, send back a null list
% of names.
if nObj == 1 && isempty(objects{1})
    objNames = [];
else
    % Every object should have a name
    objNames = cell(nObj, 1);
    for ii = 1:nObj
        if ~checkfields(objects{ii},'name')
            error('Missing object name.  %s\n',objType);
        end
        objNames{ii} = objects{ii}.name;
    end
end

end
