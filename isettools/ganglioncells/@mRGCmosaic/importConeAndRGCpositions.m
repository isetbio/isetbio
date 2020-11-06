function [coneRFpositionsMicrons, ...
          coneRFpositionsDegs, ...
          rgcRFpositionsMicrons,  ...
          rgcRFpositionsDegs, ...
          extraDegsForRGCSurround] = importConeAndRGCpositions(...
                    sourceLatticeSizeDegs, ...
                    targetEccentricityDegs, ...
                    targetSizeDegs,...
                    whichEye)
% Import positions of cones and midget RGCs from precomputed large lattices
%
% Syntax:
%   [coneRFpositionsMicrons, coneRFpositionsDegs, ...
%          rgcRFpositionsMicrons,  rgcRFpositionsDegs, ...
%          extraDegsForRGCSurround] = importConeAndRGCpositions(...
%               sourceLatticeSizeDegs, ...
%               targetEccentricityDegs, ...
%               targetSizeDegs, whichEye)
%
% Description:
%    We have precomputed large, eccentricity-varying lattices of cones and 
%    midget RGCs. These data live in mat files under: 
%           isettools/ganglioncells/data/lattices.
%    This function crops positions from the full mosaic to obtain
%    positions for mosaics of the target eccentricity and target size.
%
% Inputs:
%    sourceLatticeSizeDegs  - The size of the pre-computed lattices
%    targetEccentricityDegs - The eccentricity of the target mosaic
%    targetSizeDegs         - The size of the target mosaic 
%    whichEye               - The eye we are interested in
%
% Outputs:
%    coneRFpositionsMicrons   - The target cone positions, in microns
%    coneRFpositionsDegs      - The target cone positions, in degrees 
%    rgcRFpositionsMicrons    - The target mRGC positions, in degrees 
%    rgcRFpositionsDegs       - The target mRGC positions, in degrees 
%    extraDegsForRGCSurround  - Extra degrees to allow for cones in the
%                               surrounds of mRGCs
%
% Optional key/value pairs:
%    none
%

% History:
%    11/06/2020  NPC   Wrote it

    % Import mRGC RF positions for passed eccentricity, size, and eye
    [rgcRFpositionsMicrons,  rgcRFpositionsDegs] = retinalattice.import.finalMRGCPositions(...
        sourceLatticeSizeDegs, targetEccentricityDegs, targetSizeDegs, whichEye);

    % Determine mean surround subregion size (2 x characteristic radius) for this eccentricity
    meanEccentricityDegs = mean(sqrt(sum(rgcRFpositionsDegs.^2,2)));
    extraDegsForRGCSurround = 2.0 * ...
        RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(meanEccentricityDegs);
    
    % Import cone RF positions, for passed eccentricity, size, and eye
    [coneRFpositionsMicrons, coneRFpositionsDegs] = retinalattice.import.finalConePositions(...
        sourceLatticeSizeDegs, targetEccentricityDegs, targetSizeDegs, whichEye);                
end


