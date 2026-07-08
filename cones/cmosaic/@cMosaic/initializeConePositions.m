function initializeConePositions(obj, overlappingConeFractionForElimination)
% INITIALIZECONEPOSITIONS Imports pre-computed cone positions and
%    crops them to the desired mosaic size. 
%
% Syntax
%   obj.initializeConePositions(overlappingConeFractionForElimination)
%
% Description:
%   This method imports a large pre-computed lattice of cone positions (in
%   microns) from a source database, centers it around the mosaic's
%   eccentricity, and then crops and prepares the data for the final mosaic
%   definition. This is the default, fast way to generate a cone mosaic,
%   bypassing the slow mesh regeneration process.
%
% Inputs
%   obj  - A cMosaic object instance.
%   overlappingConeFractionForElimINATION - Scalar (0..1] or [].
%          Fraction of cone diameter used as a threshold for
%          identifying and eliminating overlapping cones. If [], no
%          overlapping cone elimination is performed. This is applied
%          to the full imported lattice before initial cropping.
%
% See also
%   cMosaic/regenerateConePositions, retinalattice.import.finalConePositions
%
% History:
%   December 2020: NPC wrote it.
%   November 2025: Gemini commented it.
%
%% 1. Import Initial Cone Positions from Pre-computed Lattice

% Import a large patch of cone positions (in microns) from the pre-computed
% lattice. The size requested is 2x the final desired size to ensure
% sufficient points for the initial cropping step. The positions are centered
% near the mosaic's 'eccentricityDegs'.
obj.coneRFpositionsMicrons = retinalattice.import.finalConePositions(...
    obj.sourceLatticeSizeDegs, obj.eccentricityDegs, obj.sizeDegs*2.0, obj.whichEye, ...
    overlappingConeFractionForElimination);

%% 2. Convert to Degrees and Perform Initial Crop (10% Buffer)

% Convert the cone positions from retinal microns to visual degrees using
% the mosaic's configured conversion function (which can be
% eccentricity-dependent).
obj.coneRFpositionsDegs = obj.distanceMicronsToDistanceDegreesForCmosaic(obj.coneRFpositionsMicrons);

% Perform an *initial* cropping to an area slightly larger than the final
% desired size (55% of half the desired size in each direction) to create
% a small buffer for the final precise crop later. This is done to speed
% up subsequent calculations.
diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
idx = find((diff(:,1) <= 0.55*obj.sizeDegs(1)) & (diff(:,2) <= 0.55*obj.sizeDegs(2)));

% Update positions to include only cones within the initial cropped region
obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);

%% 3. Compute Cone Apertures and Spacings (Used for Final Crop)

% Compute the cone apertures (light collecting areas) and spacings based
% on the current cone positions and the eccentricity-varying flags set in
% the cMosaic object. This information is needed before the final crop.
lowOpticalImageResolutionWarning = true;
obj.computeConeApertures(lowOpticalImageResolutionWarning);

%% 4. Perform Final Precise Crop
% Perform the final, precise crop of the cone positions and other related
% cone properties to exactly match the desired mosaic size ('sizeDegs')
% centered at 'eccentricityDegs'.
obj.cropMosaicDataForDesiredROI();
end

%{

Original code.

function initializeConePositions(obj, overlappingConeFractionForElimination)
% Import cone positions
% 
% Syntax
%   cMosaic.initializeConePositions(coneOverlapThreshold)
%
% Inputs
%   obj  -  A cMosaic object
%   overlappingConeFractionForElimination - Some unexplained cone overlap
%          parameter that I will try to document.
%
% , for passed eccentricity, size, and
% eye. Get a little larger region and crop after we compute the
% cone spacing below
%
%
% See also
%   cMosaic
%
% NC
%

%%
obj.coneRFpositionsMicrons = retinalattice.import.finalConePositions(...
    obj.sourceLatticeSizeDegs, obj.eccentricityDegs, obj.sizeDegs*2.0, obj.whichEye, ...
    overlappingConeFractionForElimination);

% Compute positions in degrees
obj.coneRFpositionsDegs = obj.distanceMicronsToDistanceDegreesForCmosaic(obj.coneRFpositionsMicrons);

% First crop area that is 10% wider than the desired area
diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
idx = find((diff(:,1) <= 0.55*obj.sizeDegs(1)) & (diff(:,2) <= 0.55*obj.sizeDegs(2)));
obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);

% Compute cone apertures and spacings
lowOpticalImageResolutionWarning = true;
obj.computeConeApertures(lowOpticalImageResolutionWarning);

% Crop data for desired ROI
obj.cropMosaicDataForDesiredROI();
end

%}