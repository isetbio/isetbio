function regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile, varargin)

    p = inputParser;
    p.addParameter('customRFspacingFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customMinRFspacing', [], @(x) (isempty(x) || isscalar(x)));
    p.parse(varargin{:});
    
    % Regenerate lattice whose FOV is large enough to encopass the desired size at the desired eccentricity
    fovDegs = (sqrt(sum(obj.eccentricityDegs.^2,2)) + max(obj.sizeDegs))*1.3;

    obj.coneRFpositionsMicrons = retinalattice.generatePatch(fovDegs, ...
        'cones', obj.whichEye, exportHistoryToFile, visualizeConvergence, obj.useParfor, maxIterations, ...
        'randomSeed', obj.randomSeed, ...
        'customDegsToMMsConversionFunction', obj.customDegsToMMsConversionFunction, ...
        'customRFspacingFunction', p.Results.customRFspacingFunction, ...
        'customMinRFspacing', p.Results.customMinRFspacing);
    
    % Convert to degs
    obj.coneRFpositionsDegs = obj.distanceMicronsToDistanceDegreesForCmosaic(obj.coneRFpositionsMicrons);
    
    % Compute cone apertures and spacings
    lowOpticalImageResolutionWarning = true;
    obj.computeConeApertures(lowOpticalImageResolutionWarning);
    
    % Crop data for desired ROI
    obj.cropMosaicDataForDesiredROI();
end
