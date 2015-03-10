function objNames = vcGetObjectNames(objType)
% Compile a list of object names from vcSESSION variable
%
%   objNames = vcGetObjectNames(objType)
%
% Example
%  vcGetObjectNames('oi')
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('objType'), objType = 'scene'; end

objects = vcGetObjects(objType);
nObj = length(objects);

% There may be just one, empty object.  In which case, send back a null
% list of names.
if nObj == 1 && isempty(objects{1})
    objNames = [];
else
    % Every object should have a name
    objNames = cell(nObj, 1);
    for ii=1:nObj
        if ~checkfields(objects{ii},'name')
            error('Missing object name.  %s\n',objType);
        end
        objNames{ii} = objects{ii}.name;
    end
end

end

