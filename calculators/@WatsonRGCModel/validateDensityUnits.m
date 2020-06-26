function validateDensityUnits(obj,densityUnits)
% Assert whether the passed densityUnits have valid value
%
% History:
%    2/11/20  NPC, ISETBIO Team     Wrote it.
    
    errorMessage = sprintf('densityUnits must be set to either ''%s'' or ''%s''.', obj.visualDegsDensityUnits, obj.retinalMMDensityUnits);
    assert(ismember(densityUnits, {obj.visualDegsDensityUnits, obj.retinalMMDensityUnits}), errorMessage);
end
