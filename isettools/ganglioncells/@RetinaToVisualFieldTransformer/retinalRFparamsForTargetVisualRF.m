function [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
          targetVisualRF, spatialSupportDegs, theCircularPSFData] = retinalRFparamsForTargetVisualRF(obj, ...
                  visualRFDoGparams, eccDegs, subjectEye, subjectID, subjectPupilDiameterMM, varargin)
% Estimate retinal RF params to achieve a target visual RF
%
% Syntax:
%   [retinalRFparamsStruct, weightsComputeFunctionHandle, targetVisualRF] = retinalRFparamsForTargetVisualRF(obj, ...
%          visualRFDoGparams, eccDegs, ...
%          subjectEye, subjectID, subjectPupilDiameterMM, varargin)
%
% Description:
%    Given a target visual RF and optics at a target eccentricity,
%    generate params for a retinal cone pooling RF whose visual
%    counterpart appoximates the target visual RF.
%
% Inputs:
%    visualRFDoGparams           
%    eccDegs
%    subjectEye
%    subjectID
%    subjectPupilDiameterMM
%
% Outputs:
%    retinalRFparams
%
% Optional key/value pairs:
%    'wavelengthSupport'     - Vector. Wavalength support for computing the
%                                      Vlambda weighted optics
%
% History:
%    July/August 2025  NPC  Wrote it

    p = inputParser;
    p.addParameter('wavelengthSupportForVLambdaPSF', 550 + 10*(-15:15), @isnumeric);
    p.addParameter('maxSpatialSupportDegs', 0.15, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 401, @isscalar);
    p.addParameter('psfCircularSymmetryMode', 'average', @(x)(ischar(x) && (ismember(x, {'average', 'bestResolution', 'worstResolution'}))));
    p.addParameter('deconvolutionMethod', 'Regularized', @(x)(ischar(x) && (ismember(x, {'Regularized', 'Wiener'}))));
    p.addParameter('surroundWeightBias', 0.03, @isscalar);
    p.addParameter('retinalConePoolingModel', 'GaussianCenterGaussianSurroundBased',  @(x)(ischar(x) && (ismember(x, {'GaussianCenterDoubleExponentSurroundBased', 'GaussianCenterGaussianSurroundBased'})))); 
    p.addParameter('minimizationDomain', 'visual', @(x)(ischar(x) && (ismember(x, {'retinal', 'visual'}))));
    p.parse(varargin{:});

    wavelengthSupport = p.Results.wavelengthSupportForVLambdaPSF;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    maxSpatialSupportDegs = p.Results.maxSpatialSupportDegs;
    deconvolutionMethod = p.Results.deconvolutionMethod;
    retinalConePoolingModel = p.Results.retinalConePoolingModel;
    minimizationDomain = p.Results.minimizationDomain;
    psfCircularSymmetryMode = p.Results.psfCircularSymmetryMode;
    surroundWeightBias = p.Results.surroundWeightBias;

    % Generate a @cMosaic object located at the target eccentricity and eye
    coneMosaicSize = max([0.5 2*maxSpatialSupportDegs]);
    cm = cMosaic(...
        'whichEye', subjectEye, ...
        'sizeDegs', [1 1] * coneMosaicSize, ...
        'eccentricityDegs', eccDegs, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', wavelengthSupport);

    
    % Generate a Vlambda weighted PSF and its circularly-averaged version
    % for the given subject, eye, pupil and eccentricity
    [thePSFData, theCircularPSFData] = obj.vLambdaWeightedPSFandOTF(cm, subjectID, ...
        subjectPupilDiameterMM, wavefrontSpatialSamples, maxSpatialSupportDegs, psfCircularSymmetryMode);
    

    fprintf(2,'Cone mosaic size: %2.3f degs \n', coneMosaicSize);
    fprintf(2, 'PSF size: %2.3f degs\n', (max(thePSFData.supportX)-min(thePSFData.supportX))/60);

    dxDegs = (thePSFData.supportX(2)-thePSFData.supportX(1))/60;
    x = 0:dxDegs:(2*maxSpatialSupportDegs);
    if (mod(numel(x),2) == 0)
        x(numel(x)+1) = x(end)+dxDegs;
    end
    
    rfSupportX = x - mean(x);
    rfSupportY = rfSupportX;
    spatialSupportDegs = [rfSupportX(:) rfSupportY(:)];

    [visualRFDoGparams.RcDegs, ...
     visualRFDoGparams.rotationDegs, ...
     visualRFDoGparams.centroidDegs, ...
     visualRFDoGparams.flatTopGaussianExponent, ...
     visualRFcenterConeMap, ...
     retinalRFcenterConeMap, ...
     anatomicalConeCharacteristicRadiusDegs] = estimateVisualRcFromNumberOfConesInRFcenter(...
                cm, visualRFDoGparams.conesNumPooledByTheRFcenter, ...
                theCircularPSFData, spatialSupportDegs);
        
    if (numel(visualRFDoGparams.RcDegs) == 2)

        if (visualRFDoGparams.conesNumPooledByTheRFcenter == 2)
            % If 2 cones dont let the ratio of minor/major axis be > maxRatio
            maxRxRyRatio = 2.0;
            if (visualRFDoGparams.RcDegs(1) < visualRFDoGparams.RcDegs(2))
                ratio = min([maxRxRyRatio visualRFDoGparams.RcDegs(2)/visualRFDoGparams.RcDegs(1)]);
                visualRFDoGparams.RcDegs(2) = ratio * visualRFDoGparams.RcDegs(1);
            else
                ratio = min([maxRxRyRatio visualRFDoGparams.RcDegs(1)/visualRFDoGparams.RcDegs(2)]);
                visualRFDoGparams.RcDegs(1) = ratio * visualRFDoGparams.RcDegs(2);
            end

        end

        % encode 2 radii as a complex number
        visualRFDoGparams.RcDegs = visualRFDoGparams.RcDegs(1) + 1j * visualRFDoGparams.RcDegs(2);
    end

    % The following 2 cone maps are only used when generating the summary figure
    RF2DData.visualRFcenterConeMap = visualRFcenterConeMap;
    RF2DData.retinalRFcenterConeMap = retinalRFcenterConeMap;



    % Generate the target visual RF (DoG with visualRFDoGparams and with RF center being the visual projection of the input cones)
    paramsVector = RetinaToVisualFieldTransformer.paramsStructToParamsVector(...
        visualRFDoGparams, 'GaussianCenterGaussianSurroundBased');

    RF2DData.visualRF = RetinaToVisualFieldTransformer.diffOfGaussiansRF(...
        paramsVector,  spatialSupportDegs);
    targetVisualRF = RF2DData.visualRF;

    switch (minimizationDomain)
        case 'retinal'
            % Estimate continous retinal pooling function by deconvolving
            % the target visual RF and then fitting a model to the
            % deconvolved retinal RF

            switch (deconvolutionMethod)
                case 'Regularized'
                    % Deconvolve the visual RF using the PSF data and the regularized filter
                    % algorithm, which searches for the optimal regularization parameter,
                    % alpha, within the passed regularizationAlphaRange
        
                    additiveNoisePower = 0;
                    regularizationAlphaRange = [1e-6 1];
                    RF2DData.retinalRF = deconvreg(RF2DData.visualRF,...
                        theCircularPSFData.data, additiveNoisePower, regularizationAlphaRange);
        
                case 'Wiener'
                    % Deconvolve using Wiener filter (optimal in terms of LMSE between
                    % real retinalRF and estimated retinal RF
                    RF2DData.retinalRF = deconvwnr(RF2DData.visualRF,theCircularPSFData.data);
        
                otherwise
                    error('Unknown deconvolutionMethod: ''%s''.', deconvolutionMethod);
            end
      
            switch (retinalConePoolingModel)
                case 'GaussianCenterGaussianSurroundBased'
                    % Fit DoG model to the retinal RF to extract the corresponding params
                    theFixedRetinalRcDegs = [];
                    if (visualRFDoGparams.conesNumPooledByTheRFcenter == 1)
                        theFixedRetinalRcDegs = anatomicalConeCharacteristicRadiusDegs;
                    end

                    [retinalRFparamsStruct, RF2DData.fittedRetinalRF] = ...
                        RetinaToVisualFieldTransformer.fitDoGModelToRF(...
                            spatialSupportDegs, RF2DData.retinalRF, visualRFDoGparams.rotationDegs, theFixedRetinalRcDegs);
                   
                    % Return a handle to the function that will take the
                    % retinal RF params and produce cone weights for a
                    % given cone mosaic
                    weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGmodelParameters;

                    % Compute cone indices and weights for the center & surround
                    % mechanism based on the retinalRFparamsStruct and the cone mosaic
                    pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        visualRFDoGparams.conesNumPooledByTheRFcenter, ...
                                        rfSupportX, rfSupportY, cm);
        
                    retinalRFparamsForComparisonToVisual = retinalRFparamsStruct;

                case 'GaussianCenterDoubleExponentSurroundBased'
                    theFixedRetinalRcDegs = [];
                    if (visualRFDoGparams.conesNumPooledByTheRFcenter == 1)
                        theFixedRetinalRcDegs = anatomicalConeCharacteristicRadiusDegs;
                    end

                    [retinalRFparamsStruct, RF2DData.fittedRetinalRF] = ...
                        RetinaToVisualFieldTransformer.fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(...
                            spatialSupportDegs, RF2DData.retinalRF, visualRFDoGparams.rotationDegs, theFixedRetinalRcDegs);
        
                    % Return a handle to the function that will take the
                    % retinal RF params and produce cone weights for a
                    % given cone mosaic
                    weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGDEmodelParameters;

                    % Compute cone indices and weights for the center & surround
                    % mechanism based on the retinalRFparamsStruct and the cone mosaic
                    pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        visualRFDoGparams.conesNumPooledByTheRFcenter, ...
                                        rfSupportX, rfSupportY, cm);

                    % The only parameter that can be compared to its visual counterpart is RcDegs
                    retinalRFparamsForComparisonToVisual = struct(...
                        'RcDegs', retinalRFparamsStruct.RcDegs ...
                        );

                otherwise
                    error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
            end % switch (retinalConePoolingModel)
    
            % Generate center and surround RF from the computed center/surround cone weights
            [retinalRFcenter2D, retinalRFsurround2D] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
                rfSupportX, rfSupportY, cm, pooledConeIndicesAndWeightsStruct);
        
            % And the full cone-pooling based retinal RF
            RF2DData.retinalConePoolingRF = retinalRFcenter2D - retinalRFsurround2D;
        
            % Convolve the cone-pooling based retinal RF to see if we get back the visual RF
            RF2DData.visualRFcorrespondingToRetinalConePoolingRF = conv2(RF2DData.retinalConePoolingRF, theCircularPSFData.data, 'same');

        case 'visual'
            % Here, we do not deconvolve the visualRF. Instead, we adjusting 
            % a retinal cone pooling model so as to minimize
            % the error between the visual projection of the
            % resulting retinal RF and the target visual RF

            % Constants
            modelConstants = struct();
            modelConstants.theConeMosaic = cm;
            modelConstants.thePSF = theCircularPSFData;
            modelConstants.rfSupportX = rfSupportX;
            modelConstants.rfSupportY = rfSupportY;
            modelConstants.retinalConePoolingModel = retinalConePoolingModel;
            modelConstants.conesNumPooledByTheRFcenter = visualRFDoGparams.conesNumPooledByTheRFcenter;
            modelConstants.centerConeRcDegs = anatomicalConeCharacteristicRadiusDegs;
            modelConstants.surroundWeightBias = surroundWeightBias;
            
            % Fit the visual RF by adjusting retinal cone pooling parameters
            [retinalRFparamsStruct, ...
             weightsComputeFunctionHandle, ...
             RF2DData.visualRFcorrespondingToRetinalConePoolingRF, ...
             retinalRFcenter2D, ...
             retinalRFsurround2D] = RetinaToVisualFieldTransformer.fitVisualRFByAdjustingRetinalPoolingParameters(...
                modelConstants, RF2DData.visualRF);

            RF2DData.retinalConePoolingRF = retinalRFcenter2D - retinalRFsurround2D;

            % This is the continuous retinal RF (derived by deconvolving the visual RF)
            % but since we are not deconvolving in this option (we are fitting directly the visual RF),
            % we compute it from the fitted params
            retinalRFparamsVector = RetinaToVisualFieldTransformer.paramsStructToParamsVector(retinalRFparamsStruct, retinalConePoolingModel);

            switch (retinalConePoolingModel)
                case 'GaussianCenterGaussianSurroundBased'
                    % ******* Here we must adjust the S/C integrated ratio for the fact that these params
                    % are derived from the Gaussian cone apertures. So to
                    % use them for the continuous model we must account for
                    % the difference in coverage, 
                    RF2DData.retinalRF = 0*RetinaToVisualFieldTransformer.diffOfGaussiansRF(retinalRFparamsVector, spatialSupportDegs);
                case 'GaussianCenterDoubleExponentSurroundBased'
                    % ******** Here we must adjust the S/C integrated ratio for the fact that these params
                    % are derived from the Gaussian cone apertures. So to
                    % use them for the continuous model we must account for
                    % the difference in coverage,
                    RF2DData.retinalRF = 0*RetinaToVisualFieldTransformer.diffOfGaussianCenterAndDoubleExponentSurround(retinalRFparamsVector, spatialSupportDegs);
            end

        otherwise
            error('Unknown minimizationDomain: ''%s''.', minimizationDomain);

    end

    % Finally, fit the derived visual RF (i.e., RF2DData.visualRFcorrespondingToRetinalConePoolingRF)
    % extract the achieved visual DoG params and contrast them to the desired visual DoG params)
    theFixedVisualRcDegs = [];
    
    [achievedVisualRFDoGparams, RF2DData.fittedVisualRF] = RetinaToVisualFieldTransformer.fitDoGModelToRF(...
        spatialSupportDegs, RF2DData.visualRFcorrespondingToRetinalConePoolingRF, visualRFDoGparams.rotationDegs, theFixedVisualRcDegs);

   
    switch (retinalConePoolingModel)
        case 'GaussianCenterGaussianSurroundBased'
              theFixedVisualRcDegs = [];
              retinalRFparamsForComparisonToVisual = RetinaToVisualFieldTransformer.fitDoGModelToRF(...
                    spatialSupportDegs, RF2DData.retinalRF, visualRFDoGparams.rotationDegs, theFixedVisualRcDegs);

        case 'GaussianCenterDoubleExponentSurroundBased'
              theFixedVisualRcDegs = [];
              retinalRFparamsForComparisonToVisual = RetinaToVisualFieldTransformer.fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(...
                    spatialSupportDegs, RF2DData.retinalRF, visualRFDoGparams.rotationDegs, theFixedVisualRcDegs);

              retinalRFparamsForComparisonToVisual = struct(...
                   'RcDegs', retinalRFparamsForComparisonToVisual.RcDegs ...
              );
    end

    % Generate results figure
    RetinaToVisualFieldTransformer.generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
        visualRFDoGparams, retinalRFparamsForComparisonToVisual, achievedVisualRFDoGparams, ...
        eccDegs, subjectID, maxSpatialSupportDegs, rfSupportX, rfSupportY, 1);
end




function [RcDegs, rfRotationDegs, centroidDegs, flatTopGaussianExponent, ...
    visualRFcenterConeMap, retinalRFcenterConeMap, ...
    anatomicalConeCharacteristicRadiusDegs] = estimateVisualRcFromNumberOfConesInRFcenter(cm, conesNumPooledByTheRFcenter, theCircularPSFData, spatialSupportDegs)

        % Compute the visual Rc based on the number of cones in the RF center 
        % their spacing and the PSF, all for the current eccentricity

        % Sort cones according to their distance from the mosaic center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');

        % Estimate mean anatomical cone aperture in the mosaic'c center
        sourceConesIndices = idx(1:6);
        maxConeApertureDegsInMosaicCenter = max(cm.coneApertureDiametersDegs(sourceConesIndices));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * maxConeApertureDegsInMosaicCenter;

        % Compute the retinal RF center cone map
        rfCenterPooledConeIndices = idx(1:conesNumPooledByTheRFcenter);
        meanRFCenterConePos = mean(cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:),1);
        conePosDegsRelativeToCenter = bsxfun(@minus, cm.coneRFpositionsDegs(rfCenterPooledConeIndices,:), meanRFCenterConePos);    
        retinalRFcenterConeMap = retinalRFfromPooledConeInputs(...
                anatomicalConeCharacteristicRadiusDegs, conePosDegsRelativeToCenter, spatialSupportDegs);


        % Convolve the retinal RF center cone map with the PSF
        visualRFcenterConeMap = conv2(retinalRFcenterConeMap, theCircularPSFData.data, 'same');

        % Fit a 2D Gaussian ellipsoid to the visually projected RF center cone map and extract
        % the characteristic radii of that Gaussian ellipsoid
        rfSupportX = spatialSupportDegs(:,1);
        rfSupportY = spatialSupportDegs(:,2);

        if (conesNumPooledByTheRFcenter == 2)
            % We are using a circularly summetric PSF, so force the
            % orientation to match the orientation of the 2 cones
            cone1RFpos = cm.coneRFpositionsDegs(rfCenterPooledConeIndices(1),:);
            cone2RFpos = cm.coneRFpositionsDegs(rfCenterPooledConeIndices(2),:);
            deltaY = cone2RFpos(2)-cone1RFpos(2);
            deltaX = cone2RFpos(1)-cone1RFpos(1);
            forcedOrientationDegs = -atan2d(deltaY, deltaX);
            globalSearch = true;
        else
            forcedOrientationDegs = [];
            globalSearch = false;
        end


        % Fit the visual RF center with a flat-top Gaussian ellipsoid
        [~, RcDegs, rfRotationDegs, flatTopGaussianExponent, ~, centroidDegs] = ...
            RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                rfSupportX, rfSupportY, visualRFcenterConeMap, ...
                'flatTopGaussian', true, ...
                'forcedOrientationDegs', forcedOrientationDegs, ...
                'globalSearch', globalSearch);

end

function RF2D = retinalRFfromPooledConeInputs(coneRcDegs, conePosDegs, spatialSupportDegs)
    
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));

    conesNumPooledByRFcenter = size(conePosDegs,1);
    for iCone = 1:conesNumPooledByRFcenter
  
        theConeApertureRF = exp(-((Xdegs-conePosDegs(iCone,1))/coneRcDegs).^2) .* ...
                            exp(-((Ydegs-conePosDegs(iCone,2))/coneRcDegs).^2);
        
        if (iCone == 1)
            RF2D = theConeApertureRF;
        else
            RF2D = RF2D + theConeApertureRF;
        end
    end
    RF2D = RF2D / max(RF2D(:));

end