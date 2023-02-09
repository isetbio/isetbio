classdef RTVF < handle


    % Public properties (read-only)
    properties (SetAccess = private)
       % The @cMosaic object
       coneMosaic;

       % The input optics params
       opticsParams;

       % The target visual params
       targetVisualRFDoGparams;

       % The exports directory
       exportsDirectory;

       % L- and M-cone weighted PSFs 
       % appropriate for the optical
       % position 
       spectrallyWeightedPSFData;

       % Visualization options
       visualizeSpectrallyWeightedPSFs;

       % Cone aperture blur kernel
       coneApertureBlurKernel;

       % The characteristic radius of the visual RF center as estimated by
       % fitting a 1D Gaussian line weighting function to the 1D profile
       % of the visual RF center
       visualRFcenterRcDegs;

       % The computed L-cone RF compute struct
       LconeRFcomputeStruct;

       % The computed M-cone RF compute struct
       MconeRFcomputeStruct;
    end

    % Constant properties
    properties (Constant, Hidden)
        Artal = 'Artal2012';
        Polans = 'Polans2015';
    end

    
    properties (Constant)
        validZernikeDataBases = {...
            RetinaToVisualFieldTransformer.Artal ...
            RetinaToVisualFieldTransformer.Polans ...
            };

        ArtalSubjectsNum = 54;
        PolansSubjectsNum = 10;

        defaultOpticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'positionDegs', [0 0], ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', 'right eye', ...
            'subjectRankingEye', 'right eye', ...
            'psfUpsampleFactor', 1.0, ...
            'wavefrontSpatialSamples', 401);

        % Default defaultTargetVisualRFDoGParams: no surround
        defaultTargetVisualRFDoGParams = struct(...
            'visualRFmodel', 'arbitraryShapedCenter_GaussianShapedSurround', ... 
            'retinalConePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights', ...
            'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'coneWeightsCompensateForVariationsInConeEfficiency', true,  ...
            'indicesOfConesPooledByTheRFcenter', [], ...
            'weightsOfConesPooledByTheRFcenter', [], ...
            'surroundToCenterRcRatio', 0.0, ...
            'surroundToCenterIntegratedSensitivityRatio', 0);
    end

    % Public methods (class interface)
    methods
        % Constructor
        function obj = RTVF(...
                theConeMosaic, ...
                theOpticsParams, ...
                theTargetVisualRFDoGparams, ...
                varargin) 

            p = inputParser;
            p.addParameter('exportsDirectory', '', @(x)(isempty(x)||ischar(x)));
            p.addParameter('visualizeSpectrallyWeightedPSFs', false, @islogical);
            p.addParameter('initialRetinalConePoolingParamsStruct', [], @(x)(isempty(x)||isstruct(x)));
            p.parse(varargin{:});

            % Handle optional inputs
            obj.exportsDirectory = p.Results.exportsDirectory;
            obj.visualizeSpectrallyWeightedPSFs = p.Results.visualizeSpectrallyWeightedPSFs;
            initialRetinalConePoolingParamsStruct = p.Results.initialRetinalConePoolingParamsStruct;
            
            % Handle empty opticsParams
            if (isempty(theOpticsParams))
                obj.opticsParams = RTVF.defaultOpticsParams;
            else
                obj.opticsParams = theOpticsParams;
            end

            % Handle empty cone mosaic
            if (isempty(theConeMosaic))
                obj.generateDefaultConeMosaic(theOpticsParams.positionDegs);
            else
                obj.coneMosaic = theConeMosaic;
            end
            obj.coneMosaic.visualize();

            % Handle empty argetVisualRFDoGparams
            if (isempty(theTargetVisualRFDoGparams))
                obj.targetVisualRFDoGparams = RTVF.defaultTargetVisualRFDoGParams;
            else
                obj.targetVisualRFDoGparams = theTargetVisualRFDoGparams;
            end
            
            % Compute spectrally weighted (L-cone, M-cone and L+M-cone weighted)
            %         PSFs, where the L- and M-w spectral eights are derived from the retinal
            %         spectral absorptance of actual L- and M-cones in the input cone
            %         mosaic at the desrired eccentricity (opticsParams.positionDegs), 
            %         taking into account macular pigment density at that location, as 
            %         specified by the input @cMosaic.
            %         This sets the obj.theSpectrallyWeightedPSFData parameter
            %         and updates the obj.opticsParams
            obj.spectrallyWeightedPSFs();

            % Generate the cone aperture blur kernel for the RF center pooled cones
            %          i.e., targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter
            %          This sets the obj.coneApertureKernel parameter and
            %          it is used to compute retinal cone maps for both the
            %          center and the surround subregions
            obj.computeConeApertureBlurKernel();

            % Compute the characteristic radius of the retinal RFcenter cone map
            %          as projected in visual space using the computed L+M-cone weighted PSF.
            obj.visualRFcenterRcDegs = obj.rfCenterCharacteristicRadiusInVisualSpace();


            % Retrieve the C&K Rs/Rc ratios
            [CronerKaplanRcRsRatios.temporalEccDegs, CronerKaplanRcRsRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

            % Max Rs/Rc ratio
            maxRsRcRatio =  max([ ...
                prctile(1./CronerKaplanRcRsRatios.val, 50), ...
                obj.targetVisualRFDoGparams.surroundToCenterRcRatio]);

            % Crop the PSF to this max spatial support
            obj.cropPSF(obj.visualRFcenterRcDegs * maxRsRcRatio);

            % See if we have initial retinal cone pooling params
            if (~isempty(initialRetinalConePoolingParamsStruct))
                initialLconeRetinalConePoolingParams = initialRetinalConePoolingParamsStruct.LconeRFcenter;
                initialMconeRetinalConePoolingParams = initialRetinalConePoolingParamsStruct.MconeRFcenter;
            else
                initialLconeRetinalConePoolingParams = [];
                initialMconeRetinalConePoolingParams = [];
            end

            % Action!
            obj.LconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                cMosaic.LCONE_ID, ...
                initialLconeRetinalConePoolingParams);

            obj.MconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                cMosaic.MCONE_ID, ...
                initialMconeRetinalConePoolingParams);
        end % 
        % Constructor
    end % public methods


    % Private methods
    methods (Access=private)

        % PSF computation
        spectrallyWeightedPSFs(obj, varargin);
        cropPSF(obj, maxSpatialSupportDegs);

        % Subregion computation
        computeConeApertureBlurKernel(obj);
        retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(obj, ...
            conePosDegs, coneApertureAreas, osLengthAttenuationFactors, coneWeights);

        % Rc
        visualRFcenterRcDegs = rfCenterCharacteristicRadiusInVisualSpace(obj);

        % Method to compute the RFcomputeStruct
        RFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeType, initialRetinalConePoolingParamsStruct);

    end % private methods


    % Class methods
    methods (Static)
        % Method to determine the centroid,and minor/major axis lengehs of a map
        %[theCentroid, theAxesLengths] = estimateGeometry(supportX, supportY, rfMap)

        % Method to find the rotation of theRF which results in the best
        % resolution (narrowest profile) along the x-axis
        [theRotatedRF, rotationDegs] = bestHorizontalResolutionRFmap(theRF, rotationDegs, debugRadonTransformAnalysis);

        % Method to compute the STF from the RF line weighting profile
        [oneSidedSpatialFrequencySupport, oneSidedSTF] = ...
            spatialTransferFunction(spatialSupportDegs, theRFprofile);

        % Method to fit the RF line-weighting profile with a Gaussian
        theFittedGaussianLineWeightingFunction = ...
            fitGaussianLineWeightingFunction(spatialSupport, theRFprofile, varargin)

        % Method to compute the line weighting function for a Gaussian
        theLineWeightingProfile = gaussianLineWeightingProfile(params, spatialSupport);

        % Method to visualize the params and ranges of a fitted model
        visualizeFittedModelParametersAndRanges(ax, modelParams);
    end

end

