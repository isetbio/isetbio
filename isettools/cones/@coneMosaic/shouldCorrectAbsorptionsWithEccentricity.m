function correctAbsorptions = shouldCorrectAbsorptionsWithEccentricity(obj)
% Determine whether we should correct absorptions with eccentricity
%
% Syntax:
%	correctAbsorptions = shouldCorrectAbsorptionsWithEccentricity(obj)
%
% Description:
%  Determine whether we should correct absorptions with eccentricity
%
% Inputs:
%     obj               - rect cone mosaic object
%
% Ouputs:
%     correctAbsorptions  - boolean, whether we should correct absorptions
%

% History:
%    06/16/18  NPC, ISETBIO Team    Wrote it

    correctAbsorptions = isa(obj, 'coneMosaicHex') && ...
                         (obj.eccBasedConeDensity) && ...
                         (obj.eccBasedConeQuantalEfficiency);
end

