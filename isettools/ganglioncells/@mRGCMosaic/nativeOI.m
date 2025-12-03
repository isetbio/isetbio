function [theOI, thePSF, theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio] = nativeOI(obj, varargin)

	p = inputParser;
    % Optional params
    p.addParameter('opticsModification', '', @(x)(isempty(x)||ismember(x,obj.validOpticsModifications)));
    p.addParameter('noLCA', [], @(x)(isempty(x)||islogical(x)));
    p.addParameter('customRefractionDiopters', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('wavefrontSpatialSamples', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('eccentricityDegs', obj.eccentricityDegs, @(x)(isnumeric(x)&&(numel(x)==2)));
    p.addParameter('visualizePSF', false, @islogical);
    p.addParameter('visualizedWavelengths', [], @(x)(isempty(x)||isnumeric(x)));

    p.parse(varargin{:});
    opticsModification = p.Results.opticsModification;
    eccentricityDegs = p.Results.eccentricityDegs;
    visualizePSF = p.Results.visualizePSF;
    noLCA = p.Results.noLCA;
    customRefractionDiopters = p.Results.customRefractionDiopters;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    visualizedWavelengths = p.Results.visualizedWavelengths;
    
	opticsParams = obj.rfSurroundConnectivityParams.opticsParamsStruct;

	if (~isempty(opticsModification))
		opticsParams.modification = opticsModification;
	end

	if (~isempty(noLCA))
		opticsParams.noLCA = noLCA;
	end

	opticsParams.refractiveErrorDiopters = 0.0;
	if (~isempty(customRefractionDiopters))
		opticsParams.refractiveErrorDiopters = customRefractionDiopters;
	end

	eccentricityRadiusDegs = sqrt(sum(eccentricityDegs.^2,2));
	mosaicMaxRadiusDegs = sqrt(sum((abs(obj.eccentricityDegs) + obj.sizeDegs).^2,2));
	if (eccentricityRadiusDegs > mosaicMaxRadiusDegs)
		error('Generated optics are at an eccentricity (%2.1f) that is outside of the max eccentricity of the mosaic (%2.1f)', ...
			eccentricityRadiusDegs, mosaicMaxRadiusDegs);
	end


	[theOI, thePSF, theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio] = RGCMosaicConstructor.helper.optics.generate(...
    	obj.inputConeMosaic, ...
    	eccentricityDegs, opticsParams, ...
        'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
    	'visualizePSF', visualizePSF, ...
    	'visualizedWavelengths',visualizedWavelengths ,...
    	'visualizeStrehlRatioOptimization', false);
end
