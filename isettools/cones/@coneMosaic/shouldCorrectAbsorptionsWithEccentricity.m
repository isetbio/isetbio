function correctAbsorptions = shouldCorrectAbsorptionsWithEccentricity(obj)
% Determine whether we should correct absorptions with eccentricity
%
% Syntax:
%	correctAbsorptions = shouldCorrectAbsorptionsWithEccentricity(obj)
%
% Description:
%  Determine whether we should correct absorptions with eccentricity. This
%  only happens when the following 3 conditions are all true:
%    - the mosaic is hexagonal 
%    - the eccBasedConeDensity flag is set to true 
%    - the eccBasedConeQuantalEfficiency flag is set to true
%
% Inputs:
%     obj                 -  cone mosaic object
%
% Ouputs:
%     correctAbsorptions  - boolean, if true, then the computeSingleFrame()
%                           method corrects the absorptions to account for
%                           changes in cone efficiency with eccentricity
%

% History:
%    06/16/18  NPC, ISETBIO Team    Wrote it

    correctAbsorptions = isa(obj, 'coneMosaicHex') && ...
                         (obj.eccBasedConeDensity) && ...
                         (obj.eccBasedConeQuantalEfficiency);
end

