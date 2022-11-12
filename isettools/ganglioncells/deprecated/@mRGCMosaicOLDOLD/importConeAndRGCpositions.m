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

    % Determine surround subregion size (2 x characteristic radius) for the max eccentricity
    maxEccentricityDegs = max(sqrt(sum(rgcRFpositionsDegs.^2,2)));
    extraDegsForRGCSurround = 2.0 * ...
        RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);
    
    % Import cone RF positions, for passed eccentricity, size, and eye
    [coneRFpositionsMicrons, coneRFpositionsDegs] = retinalattice.import.finalConePositions(...
        sourceLatticeSizeDegs, targetEccentricityDegs, targetSizeDegs, whichEye);   
    
    fprintf('Imported %d cone and %d RGC positions covering an area of %2.1f x %2.1f degs at ecc = (%2.1f,%2.1f) degs\n', ...
         size(coneRFpositionsMicrons,1), size(rgcRFpositionsMicrons,1), ...
         targetSizeDegs(1), targetSizeDegs(2), targetEccentricityDegs(1), targetEccentricityDegs(2));

end


