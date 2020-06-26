function validateEccUnits(obj,eccUnits)
% Assert whether the passed eccUnits have valid value
%
% History:
%    2/11/20  NPC, ISETBIO Team     Wrote it.
    
    errorMessage = sprintf('eccUnits must be set to either ''%s'' or ''%s''.', obj.visualDegsEccUnits, obj.retinalMMEccUnits);
    assert(ismember(eccUnits, {obj.visualDegsEccUnits, obj.retinalMMEccUnits}), errorMessage);
end
