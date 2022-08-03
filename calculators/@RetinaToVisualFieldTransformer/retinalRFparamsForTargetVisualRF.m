function [retinalRFparamsStruct, weightsComputeFunctionHandle, targetVisualRF, theCircularPSFData] = retinalRFparamsForTargetVisualRF(obj, ...
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
    p.addParameter('psfCircularSymmetryMode', 'average', @(x)(ischar(x) && (ismember(x, {'average', 'best', 'worse'}))));
    p.addParameter('deconvolutionMethod', 'Regularized', @(x)(ischar(x) && (ismember(x, {'Regularized', 'Wiener'}))));
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

    % Generate a @cMosaic object located at the target eccentricity and eye
    fprintf(2,'Cone mosaic size: %2.3f degs \n', max([0.5 maxSpatialSupportDegs]));
    pause(1);
    cm = cMosaic(...
        'whichEye', subjectEye, ...
        'sizeDegs', [1 1] * max([0.5 maxSpatialSupportDegs]), ...
        'eccentricityDegs', eccDegs, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', wavelengthSupport);

    
    % Generate a Vlambda weighted PSF and its circularly-averaged version
    % for the given subject, eye, pupil and eccentricity
    [thePSFData, theCircularPSFData] = obj.vLambdaWeightedPSFandOTF(cm, subjectID, ...
        subjectPupilDiameterMM, wavefrontSpatialSamples, maxSpatialSupportDegs, psfCircularSymmetryMode);
    

    % If an RcDegs is not specified, compute it from the conesNumPooledByTheRFcenter
    % and the optics
    if (isempty(visualRFDoGparams.RcDegs))
        [visualRFDoGparams.RcDegs, visualRFcenterConeMap, ...
         retinalRFcenterConeMap, anatomicalConeCharacteristicRadiusDegs] = ...
            RetinaToVisualFieldTransformer.estimateVisualRcFromNumberOfConesInRFcenter(...
                cm, visualRFDoGparams.conesNumPooledByTheRFcenter, theCircularPSFData);
        
        % The following 2 cone maps are only used when generating the summary figure
        RF2DData.visualRFcenterConeMap = visualRFcenterConeMap;
        RF2DData.retinalRFcenterConeMap = retinalRFcenterConeMap;
    end


    % Generate the target visual RF (DoG with visualRFDoGparams and with RF center being the visual projection of the input cones)
    spatialSupportDegs = [thePSFData.supportX(:) thePSFData.supportY(:)]/60;
    paramsVector = RetinaToVisualFieldTransformer.paramsStructToParamsVector(visualRFDoGparams, retinalConePoolingModel);

    RF2DData.visualRF = RetinaToVisualFieldTransformer.diffOfGaussiansRF(paramsVector,  spatialSupportDegs);
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
                    [RF2DData.retinalRF, lagra] = deconvreg(RF2DData.visualRF,...
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
                        RetinaToVisualFieldTransformer.fitDoGModelToRF(spatialSupportDegs, RF2DData.retinalRF, theFixedRetinalRcDegs);
                   
                    % Return a handle to the function that will take the
                    % retinal RF params and produce cone weights for a
                    % given cone mosaic
                    weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGmodelParameters;

                    % Compute cone indices and weights for the center & surround
                    % mechanism based on the retinalRFparamsStruct and the cone mosaic
                    pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        visualRFDoGparams.conesNumPooledByTheRFcenter, ...
                                        thePSFData.supportX, thePSFData.supportY, cm);
        
                    retinalRFparamsForComparisonToVisual = retinalRFparamsStruct;

                case 'GaussianCenterDoubleExponentSurroundBased'
                    theFixedRetinalRcDegs = [];
                    if (visualRFDoGparams.conesNumPooledByTheRFcenter == 1)
                        theFixedRetinalRcDegs = anatomicalConeCharacteristicRadiusDegs;
                    end

                    [retinalRFparamsStruct, RF2DData.fittedRetinalRF] = ...
                        RetinaToVisualFieldTransformer.fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(...
                            spatialSupportDegs, RF2DData.retinalRF, theFixedRetinalRcDegs);
        
                    % Return a handle to the function that will take the
                    % retinal RF params and produce cone weights for a
                    % given cone mosaic
                    weightsComputeFunctionHandle = @RetinaToVisualFieldTransformer.retinalConeWeightsFromDoGDEmodelParameters;

                    % Compute cone indices and weights for the center & surround
                    % mechanism based on the retinalRFparamsStruct and the cone mosaic
                    pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
                                        retinalRFparamsStruct, ...
                                        visualRFDoGparams.conesNumPooledByTheRFcenter, ...
                                        thePSFData.supportX, thePSFData.supportY, cm);

                    % The only parameter that can be compared to its visual counterpart is RcDegs
                    retinalRFparamsForComparisonToVisual = struct(...
                        'RcDegs', retinalRFparamsStruct.RcDegs ...
                        );

                otherwise
                    error('Unknown retinalConePoolingModel: ''%s''.', retinalConePoolingModel);
            end % switch (retinalConePoolingModel)
    
            % Generate center and surround RF from the computed center/surround cone weights
            [retinalRFcenter2D, retinalRFsurround2D] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
                thePSFData.supportX, thePSFData.supportY, cm, pooledConeIndicesAndWeightsStruct);
        
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
            modelConstants.retinalConePoolingModel = retinalConePoolingModel;
            modelConstants.conesNumPooledByTheRFcenter = visualRFDoGparams.conesNumPooledByTheRFcenter;
            modelConstants.centerConeRcDegs = anatomicalConeCharacteristicRadiusDegs;
           
            % Fit the visual RF by adjusting retinal cone pooling parameters
            [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
             RF2DData.visualRFcorrespondingToRetinalConePoolingRF, ...
             retinalRFcenter2D, retinalRFsurround2D] = RetinaToVisualFieldTransformer.fitVisualRFByAdjustingRetinalPoolingParameters(...
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
                    RF2DData.retinalRF = RetinaToVisualFieldTransformer.diffOfGaussiansRF(retinalRFparamsVector, spatialSupportDegs);
                case 'GaussianCenterDoubleExponentSurroundBased'
                    % ******** Here we must adjust the S/C integrated ratio for the fact that these params
                    % are derived from the Gaussian cone apertures. So to
                    % use them for the continuous model we must account for
                    % the difference in coverage,
                    RF2DData.retinalRF = RetinaToVisualFieldTransformer.diffOfGaussianCenterAndDoubleExponentSurround(retinalRFparamsVector, spatialSupportDegs);
            end

        otherwise
            error('Unknown minimizationDomain: ''%s''.', minimizationDomain);

    end

    % Finally, fit the derived visual RF (i.e., RF2DData.visualRFcorrespondingToRetinalConePoolingRF)
    % extract the achieved visual DoG params and contrast them to the desired visual DoG params)
    theFixedVisualRcDegs = [];
    [achievedVisualRFDoGparams, RF2DData.fittedVisualRF] = RetinaToVisualFieldTransformer.fitDoGModelToRF(...
        spatialSupportDegs, RF2DData.visualRFcorrespondingToRetinalConePoolingRF, theFixedVisualRcDegs);

   
    switch (retinalConePoolingModel)
        case 'GaussianCenterGaussianSurroundBased'
              theFixedVisualRcDegs = [];
              retinalRFparamsForComparisonToVisual = RetinaToVisualFieldTransformer.fitDoGModelToRF(...
                    spatialSupportDegs, RF2DData.retinalRF, theFixedVisualRcDegs);

        case 'GaussianCenterDoubleExponentSurroundBased'
              theFixedVisualRcDegs = [];
              retinalRFparamsForComparisonToVisual = RetinaToVisualFieldTransformer.fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(...
                    spatialSupportDegs, RF2DData.retinalRF, theFixedVisualRcDegs);

              retinalRFparamsForComparisonToVisual = struct(...
                   'RcDegs', retinalRFparamsForComparisonToVisual.RcDegs ...
              );
    end

    % Generate results figure
    RetinaToVisualFieldTransformer.generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
        visualRFDoGparams, retinalRFparamsForComparisonToVisual, achievedVisualRFDoGparams, ...
        eccDegs, subjectID, maxSpatialSupportDegs, 1);
end
