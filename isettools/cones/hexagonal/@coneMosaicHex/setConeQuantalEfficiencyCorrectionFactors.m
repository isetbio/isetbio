function setConeQuantalEfficiencyCorrectionFactors(obj, correctionFactors)
% Set the cone quantal efficiency correction factors
%
% Syntax:
%    setConeQuantalEfficiencyCorrectionFactors(obj, correctionFactors)
%
% Description:
%    Called by the static method coneMosaicHex.computeConeEfficiencyCorrectionFactors()
%
% Inputs:
%    obj - The cone mosaic hex object
%    correctionFactors - The correction factors
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    06/15/2018  NPC  ISETBIO TEAM . Wrote it
%
    obj.coneEfficiencyCorrectionFactors =  correctionFactors;
end