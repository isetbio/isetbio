classdef CronerKaplanRGCModel < handle
    % Create a CronerKaplan RGC Model
    
    % References:
    %    Croner&Kaplan (1994). 'Receptive fields of P and M Ganglion cells 
    %    across the primate retina.',Vis. Res, (35)(1), pp.7-24
    % History:
    %    11/8/19  NPC, ISETBIO Team     Wrote it.
    
    
    properties (SetAccess = private)
        % Digitized data from Figure 4 & 5
        centerData;
        surroundData;
        
        % Digitized data from Figures 4,6 & 11
        surroundCenterRatiosData;
        
        % Digitized data form Figure 3
        recordedEccentricitiesData;
        
        % Synthesized data
        synthesizedData;
        
        % Model of center radius with eccentricity
        centerRadiusFunction;
        centerRadiusThreshold;
        centerRadiusParams;
        centerRadiusParamsSE;
        
        % Model of surround radius with eccentricity
        surroundRadiusFunction;
        surroundRadiusThreshold;
        surroundRadiusParams;
        surroundRadiusParamsSE;
        
        % Model of center sensitivity with center radius
        centerPeakSensitivityFunction;
        centerPeakSensitivityParams;
        centerPeakSensitivityParamsSE;
        
        % Model of surround sensitivity with surround radius
        surroundPeakSensitivityFunction;
        surroundPeakSensitivityParams;
        surroundPeakSensitivityParamsSE;
        
        % Valid quadrant names for the Polans wavefront-based optics
        validPolansQuadrants = {'horizontal', 'superior', 'inferior'};
        validPolansSubjectIDs = 1:10;
        
        % Directory with psf deconvolution results
        psfDeconvolutionDir;
        
        % Eccs for which we conduced deconvolution 
        % (see obj.generateDeconvolutionFiles())
        deconvolutionEccs;
        
        synthesisOptions;
        
        plotlabOBJ;
    end
    
    properties (Constant)
        % Eccentricities for which we have generated deconvolution data files
        % (via generateDeconvolutionFilesForMidgetRGCs() )
        defaultDeconvolutionEccs = [0 0.1 0.2 0.3 0.5 0.7 1 1.5 2 2.5 3:21];
        
        % Paper formulas
        % Surround radius from ecc (Figure 4 caption)
        surroundRadiusFromEccDegs = @(eccDegs) (0.203 * abs(eccDegs).^0.472);
        
        % Surround radius from center radius
        surroundRadiusFromCenterRadiusDegs = @(centerRadiusDegs) (centerRadiusDegs * 6.7);
        
        % Surround peak sensitivity from surround radius (Figure 5 caption)
        surroundPeakSensitivityFromSurroundRadiusDegs = @(surroundRadiusDegs) (0.128 * surroundRadiusDegs.^(-2.147));
        
        % Center peak sensitivity from center radius (Figure 5 caption)
        centerPeakSensitivityFromCenterRadiusDegs = @(centerRadiusDegs) (0.391 * centerRadiusDegs.^(-1.850));
        
        % Surround:Center integrated sensitivity ratio (Figure 11 caption)
        surroundToCenterIntegratedVisualSensitivityRatiosFromEccDegs = @(eccDegs)  (0.466 + 0.007 * abs(eccDegs));
    end
    
    methods
        % Constructor
        function obj = CronerKaplanRGCModel(varargin) 

            % Parse input
            p = inputParser;
            p.addParameter('generateAllFigures', true, @islogical);
            p.addParameter('instantiatePlotLab', true, @islogical);
            p.addParameter('deconvolutionEccs', CronerKaplanRGCModel.defaultDeconvolutionEccs, @isnumeric);
            p.addParameter('dataSetToFit', 'medians', @(x)(ismember(x, {'medians', 'raw', 'paperFormulas'})));
            p.parse(varargin{:});
            
            obj.psfDeconvolutionDir = strrep(fileparts(which(mfilename())), ...
                '@CronerKaplanRGCModel', 'VisualToRetinalCorrectionData/DeconvolutionData');
            
            obj.deconvolutionEccs = p.Results.deconvolutionEccs;
            
            obj.loadRawData();
            obj.fitModel('dataset', p.Results.dataSetToFit);
            
            obj.synthesisOptions = struct( ...
                'randomizeCenterRadii', true, ...
                'randomizeCenterSensitivities', true, ...
                'randomizeSurroundRadii', true, ...
                'randomizeSurroundSensitivities', true);
            
            if (p.Results.instantiatePlotLab)
                obj.setupPlotLab(0, 14, 14);
            end
            
            if (p.Results.generateAllFigures)
                obj.plotDigitizedData();
            end
        end
        
        % Fit the model to a data set, either 'medians', or 'raw'
        fitModel(obj, varargin);
        
        % Method to synthesize data for a sample of eccentricity values
        synthesizeData(obj, eccDegs, synthesisOptions);
        
        % Method to plot different aspects of the synthesized data
        [hFig1, hFig2, hFig3, hFig4] = plotSynthesizedData(obj);
        
        % Method to simulate the Croner&Kaplan results
        simulateCronerKaplanResults(obj, varargin);
        
        % Assemble deconvolution file name
        dataFileName = deconvolutionDataFileName(obj, patchEccRadiusDegs, ...
            imposedRefractionErrorDiopters, subregionName);
        
        % Generate the Gaussian-PSF deconvolution analysis data files
        generateDeconvolutionFiles(obj, deconvolutionOpticsParams, subregion, varargin);

        % Perform the Gaussian-PSF deconvolution analysis for the RF center
        deconvolutionStruct = performDeconvolutionAnalysisForRFcenter(obj, conesNumInRFcenterExamined, ...
            sensitivityRangeOverWhichToMatchSFtuning, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
            quadrantName, subjectID, patchEccRadiusDegs);
        
        % Perform the Gaussian-PSF deconvolution analysis for the RF surround
        deconvolutionStruct = performDeconvolutionAnalysisForRFsurround(obj, deconvolutionStructForRFcenter, ...
            examinedConesNumInRFCenter, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
            quadrantName, subjectID, patchEccRadiusDegs);
        
        % Generate the deconvolution model (operates on the output of
        % performGaussianConvolutionWithPolansPSFanalysis()) - no printing
        deconvolutionModel = computeDeconvolutionModel(obj, deconvolutionOpticsParams);
        
        % Method to generate retinal RF params given the retinal center radius
        % and eccentricity as inputs. This uses (via computeDeconvolutionModel()),
        % the Gaussian-PSF convolution data generated by performGaussianConvolutionWithPolansPSFanalysis().
        % It is to be used with mRGC mosaics whose centers are determined by connectivity to an underlying
        % cone mosaic.
        synthesizedRFParams = synthesizeRetinalRFparamsConsistentWithVisualRFparams(obj, ...
            retinalCenterInputConesNum, retinalCenterMicrons, deconvolutionOpticsParams);
    end
    
    methods (Static)
        plotSensitivities(theAxes, d, model, pointSize, color,displayYLabel, theLabel);
        
        plotRadii(theAxes, d, model, pointSize, color, displayYLabel, theLabel);
        
        
        % Generate a cone mosaic for performing the deconvolution analysis
        [theConeMosaic, theConeMosaicMetaData] = generateConeMosaicForDeconvolution(patchEcc, patchSize, varargin);
        
        % Generate Polans optics for performing the deconvolution analysis
        [theOptics, theFullPSF, thePSFsupportDegs] = generatePolansOpticsForDeconcolution(PolansSubjectID, ...
            imposedRefractionErrorDiopters, pupilDiameterMM , wavelengthSampling, micronsPerDegree, patchEcc, varargin);
        
        % Generate the cone aperture profile for the deconvolution analysis
        [coneApertureProfile,  theConeApertureSupportDegs] = generateConeApertureProfileForDeconvolution(thePSFsupportDegs, coneAperturesDegs);
        
        % Interpolate and zero-pad the PSF
        [thePSFHR, thePSFsupportDegsHR] = interpolatePSF(thePSF, thePSFsupportDegs, upsampleFactor, paddingMarginDegs);
        
        % Integrate retinal/visual cone image within the apertures of cones for the deconvolution analysis
        [retinalConeActivations, visualConeActivations, ...
         withinConeAperturesRetinalConeImage, ...
         withinConeAperturesVisualConeImage, ...
         withinConeAperturesGrid] = integrateConeImageWithinConeAperturesForDeconvolution(...
                    retinalConeImage, visualConeImage, coneMask, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

        % Fit retinal/visual activation map for the deconvolution analysis
        [rfSigmas, rfGain, ...
         hiResConeActivationMap, hiResPSFsupportDegs] = fitActivationMapForDeconvolution(...
                   functionName, coneActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

        % Fit an elliptical 2D Gaussian to a RF
        [fittedParams,  rfFunction] = fitElliptical2DGausianToRF(functionName, X, Y, RF, deltaX, minSigma, center);
  
        % Visualize the retinal and visual RF subregion fits for the deconvolution analysis
        visualizeFitsForDeconvolution(inputConeCombinationIndex, thePSFsupportDegs, thePSF, ...
                conePosDegs, retinalConeImage, visualConeImage, ...
                retinalConeActivations, visualConeActivations, ...
                hiResPSFsupportDegs, hiResRetinalConeActivationMap, hiResVisualConeActivationMap, ...
                rfSigmasRetinal, rfGainRetinal, rfSigmasVisual, rfGainVisual, ...
                visualizedSpatialSupportRange);
            
        % Generate Polans PSF at desired eccentricitiy
        [hEcc, vEcc, thePSFs, thePSFsupportDegs, theOIs] = psfsAtEccentricity(goodSubjects, ...
            imposedRefractionErrorDiopters, desiredPupilDiamMM, wavelengthsListToCompute, ...
            micronsPerDegree, wavefrontSpatialSamples, eccXrange, eccYrange, deltaEcc, varargin);
        
        data = quadrantData(allQuadrantData, quadrantsToAverage, quadrantsComputed, subjectsToAverage, subjectsComputed);
        
        % Print the info contained in a deconvolutionStruct
        printDeconvolutionStruct(deconvolutionStruct);
        
        plotDeconvolutionModel(deconvolutionModel);
    end
    
    methods (Access=private)
        setupPlotLab(obj, mode, figWidthInches, figHeightInches);
        
        % Method to validate the deconvolutionOpticsParams
        validateDeconvolutionOpticsParams(obj,deconvolutionOpticsParams);
    end
    
end

