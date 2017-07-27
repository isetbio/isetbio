function obj = gatherStruct(obj)
% Gather distributed struct to current working directory
%   function obj = gatherStruct(obj)
%
% Input:
%   obj  - variable or structure, for structure, we will gather recursively
%          for all its sub-field
%
% Output:
%   obj  - gathered obj
%
% Example:
%   scene = sceneCreate;
%   scene = gather(scene);
%
% Note:
%   This function is only useful in context of gpu computing. If there's
%   no distributed field (e.g. gpuArray) in obj, the output will be the
%   same as input
%
% See also:
%   gather, vcAddObject, vcAddAndSelectObject
%
% (HJ) ISETBIO TEAM, 2014

if notDefined('obj'), obj = []; return; end

if isstruct(obj)
    fNames = fieldnames(obj);
    for ii = 1:length(fNames)
        obj.(fNames{ii}) = gatherStruct(obj.(fNames{ii}));
    end
elseif isa(obj, 'gpuArray')
    obj = gather(obj);
end

end