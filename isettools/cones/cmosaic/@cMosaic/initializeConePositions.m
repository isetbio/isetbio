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

