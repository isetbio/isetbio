function [theCenterConnectedMRGCMosaicFullFileName, theIntermediateConnectivityStageMetaDataFile, ...
          theCenterConnectedMRGCMosaicFileName] = centerConnectedMosaic(...
            whichEye, eccentricityDegs, sizeDegs, rfCenterConnectivityParams, sourceLatticeSizeDegs,  varargin)
    
    % Parse input
    p = inputParser;

    % Required params
    p.addRequired('whichEye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
    p.addRequired('eccentricityDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
    p.addRequired('sizeDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
    p.addRequired('sourceLatticeSizeDegs', @isscalar);
    p.addRequired('rfCenterConnectivityParams', @(x)rfCenterConnectivityValidationFunction(x));

    % Optional params
    % Custom Degs <-> MM conversion functions
    p.addParameter('customLMSconeDensities', [], @(x)(isempty(x)||((isnumeric(x))&&(numel(x)==3))));
    p.addParameter('customDegsToMMsConversionFunction', @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x), @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customMMsToDegsConversionFunction', @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x), @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('randomSeed', [], @(x)(isempty(x) || isscalar(x)));

    p.parse(whichEye, eccentricityDegs, sizeDegs, ...
            sourceLatticeSizeDegs, rfCenterConnectivityParams, ...
            varargin{:});

    % Generate the filename where the center-connected RGC mosaic will be saved
    pCenter = p.Results;
    pCenter.spatialChromaticUniformityTradeoff = rfCenterConnectivityParams.spatialChromaticUniformityTradeoff;
    [theCenterConnectedMRGCMosaicFullFileName, ...
     theIntermediateConnectivityStageMetaDataFile, ...
     theCenterConnectedMRGCMosaicFileName] = RGCMosaicConstructor.filepathFor.exportedMosaicFileName(pCenter, 'center connected');


    % STEP 1. Generate the input cone mosaic
    inputConeMosaicParams = p.Results;
    [theInputConeMosaic, extraDegsToSupportRGCSurrounds] = RGCMosaicConstructor.compute.inputConeMosaic(inputConeMosaicParams);
    theInputConeMosaic.visualize();


    % STEP 2. Generate the RGC RF centers and the cones-RF centers connectivity matrix
    theRGCMosaicComponents = RGCMosaicConstructor.compute.componentsForRFcenterConnectivity(...
        theInputConeMosaic, extraDegsToSupportRGCSurrounds, ...
        rfCenterConnectivityParams, theIntermediateConnectivityStageMetaDataFile);

    % Add the mosaic params
    theRGCMosaicComponents.whichEye = whichEye;
    theRGCMosaicComponents.eccentricityDegs = eccentricityDegs;
    theRGCMosaicComponents.sizeDegs = sizeDegs;

    % Generate mRGC mosaic with the center-connection components
    theMRGCMosaic = mRGCMosaic('withComponents', theRGCMosaicComponents);

    % Generate the filename for saving the center-connected RGC mosaic
    fprintf('Saving center-connected MRGCMosaic to %s\n',theCenterConnectedMRGCMosaicFullFileName);
    save(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic', '-v7.3');
end


function s = rfCenterConnectivityValidationFunction(x)
    c(1) = isstruct(x);
    c(numel(c)+1) = isfield(x, 'optimizationCenter');
    c(numel(c)+1) = isfield(x, 'spatialChromaticUniformityTradeoff');
    c(numel(c)+1) = isfield(x, 'coneTypesToBeConnected');
    c(numel(c)+1) = isfield(x, 'maxNeighborNormDistance');
    c(numel(c)+1) = isfield(x, 'maxNeighborsNum');
    c(numel(c)+1) = isfield(x, 'localSpacingFromCurrentCentroids');
    c(numel(c)+1) = isfield(x, 'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs');
    c(numel(c)+1) = isfield(x, 'visualizeIntermediateStagesOfCenterConnectivityOptimization');
    c(numel(c)+1) = isfield(x, 'saveIntermediateStagesOfCenterConnectivityOptimization');

    if (~all(c))
        x
        error('Invalid rfCenterConnectivityStruct');
    else
       s = true;
    end
end
