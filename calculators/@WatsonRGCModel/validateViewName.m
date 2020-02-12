function validateViewName(obj,viewName)
% Assert whether the passed view has valid name
%
% History:
%    2/11/20  NPC, ISETBIO Team     Wrote it.
    
    errorMessage = sprintf('viewName must be set to either ''%s'', ''%s'' or ''%s''.', obj.rightEyeVisualField, obj.rightEyeRetina, obj.leftEyeRetina);
    assert(ismember(viewName, {obj.rightEyeVisualField, obj.rightEyeRetina, obj.leftEyeRetina}), errorMessage);
end

        
