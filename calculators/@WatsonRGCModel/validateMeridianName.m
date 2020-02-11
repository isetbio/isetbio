function validateMeridianName(obj,meridianName)
% Assert whether the passed meridianName has valid value
%
% History:
%    2/11/20  NPC, ISETBIO Team     Wrote it.
    
    errorMessage = sprintf('meridian name must be set to either ''%s'', ''%s'', ''%s'', or ''%s''.', ...
        obj.enumeratedMeridianNames{1}, obj.enumeratedMeridianNames{2}, ...
        obj.enumeratedMeridianNames{3}, obj.enumeratedMeridianNames{4});
    assert(ismember(meridianName, obj.enumeratedMeridianNames), errorMessage);
end
