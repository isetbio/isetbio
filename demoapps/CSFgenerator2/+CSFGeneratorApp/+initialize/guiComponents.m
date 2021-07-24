function guiComponents(app)
    % Initialize colors for different aspects of the app
    initializeColors(app);
    
    % Initialize the status fields
    initializeStatusFields(app);

    % Initialize the region of interest GUI components
    initializeROIGUIComponents(app);
    
    % Initialize the stimulus GUI components
    initializeStimulusGUIComponents(app);
    
    % Initialize the optics GUI components
    initializeOpticsGUIComponents(app);
    
    % Initialize the cone mosaic GUI components
    initializeConeMosaicGUIComponents(app);
    
    % Initialize the computational observer GUI components
    initializeComputationalObserverGUIComponents(app);
end


function initializeColors(app)
    app.colors = containers.Map();
    app.colors('good message background') = [0.2 0.3 0.7];
    app.colors('good message foreground') = [1 1 1];
    app.colors('problem message background') = [1 0 0];
    app.colors('problem message foreground') = [1 1 1];
end


function initializeStatusFields(app)

    s = struct(...
        'text', ' ', ...
        'fontColor', app.colors('good message foreground'), ...
        'backgroundColor', app.colors('good message background'), ...
        'fontWeight', 'normal');
    
    app.statusMessages = containers.Map();
    app.statusMessages('region of interest (visual field)') = s;
    app.statusMessages('stimulus') = s;
    app.statusMessages('tasks') = s;
    app.statusMessages('performance assessment') = s;
    app.statusMessages('optics') = s;
    app.statusMessages('cone mosaic') = s;
    app.statusMessages('fixational eye movements') = s;
    app.statusMessages('mRGC mosaic') = s;
    app.statusMessages('pRGC mosaic') = s;
    app.statusMessages('computational observer') = s;
    
    CSFGeneratorApp.render.statusField(app,'A', 'region of interest (visual field)');
    CSFGeneratorApp.render.statusField(app,'B', 'optics'); 
end

function initializeStimulusGUIComponents(app)
    % Stimulus presentation mode, duration, and temporal frequency
    CSFGeneratorApp.decode.stimulusPresentationModeDropDown(app, 'valueToSlider', app.stimParams.presentationMode);
    CSFGeneratorApp.decode.stimulusDurationSpinner(app, 'valueToSlider', app.stimParams.durationSec);
    CSFGeneratorApp.decode.stimulusTemporalFrequencySpinner(app, 'valueToSlider', app.stimParams.temporalFrequencyHz);
    
    % Stimulus wavelength
    CSFGeneratorApp.decode.stimulusWavelengthSupportStepSizeSpinner(app, 'valueToSlider', app.stimParams.wavelengthSupportStepSize);
    CSFGeneratorApp.decode.stimulusWavelengthSupportMaxSpinner(app, 'valueToSlider', app.stimParams.wavelengthSupportMax);
    CSFGeneratorApp.decode.stimulusWavelengthSupportMinSpinner(app, 'valueToSlider', app.stimParams.wavelengthSupportMin);
    
    % Stimulus LMS contrast
    CSFGeneratorApp.decode.stimulusLconeContrastSpinner(app, 'valueToSlider', app.stimParams.LconeContrast);
    CSFGeneratorApp.decode.stimulusMconeContrastSpinner(app, 'valueToSlider', app.stimParams.MconeContrast);
    CSFGeneratorApp.decode.stimulusSconeContrastSpinner(app, 'valueToSlider', app.stimParams.SconeContrast);
    
    % Stimulus mean luminance
    CSFGeneratorApp.decode.stimulusMeanLuminanceSpinner(app, 'valueToSlider', app.stimParams.meanLuminanceCdM2);
    
    % Stimulus position
    CSFGeneratorApp.decode.stimulusSpatialPositionXSpinner(app, 'valueToSlider', app.stimParams.positionDegs(1));
    CSFGeneratorApp.decode.stimulusSpatialPositionYSpinner(app, 'valueToSlider', app.stimParams.positionDegs(2));
    CSFGeneratorApp.decode.stimulusMosaicCenteredCheckBox(app, 'valueToSlider', app.stimParams.mosaicCenteredPosition);
    
    % Stimulus resolution and size
    CSFGeneratorApp.decode.stimulusResolutionSpinner(app, 'valueToSlider', app.stimParams.resolutionPixels);
    CSFGeneratorApp.decode.stimulusSizeSpinner(app, 'valueToSlider', app.stimParams.sizeDegs);
    
    % Stimulus spatial envelope, orientation, spatial frequency and spatial phase
    CSFGeneratorApp.decode.stimulusSpatialEnvelopeDropDown(app, 'valueToSlider', app.stimParams.spatialEnvelope);
    CSFGeneratorApp.decode.stimulusSpatialFrequencySpinner(app, 'valueToSlider', app.stimParams.spatialFrequencyCPD);
    CSFGeneratorApp.decode.stimulusSpatialPhaseSpinner(app, 'valueToSlider', app.stimParams.spatialPhaseDegs);
    CSFGeneratorApp.decode.stimulusOrientationSpinner(app, 'valueToSlider', app.stimParams.orientationDegs);
end

function initializeOpticsGUIComponents(app)
    % The visualized wavelength slider
    app.opticsVisualizedWavelengthSlider.Limits = [400 700];
    app.opticsVisualizedWavelengthSlider.MajorTicks = 400:50:700;
    app.opticsVisualizedWavelengthSlider.MinorTicks = 400:10:700;
    CSFGeneratorApp.decode.opticsWavelengthSlider(app, 'valueToSlider', app.opticsParams.visualizedWavelength);
    
    % The subject dataset
    CSFGeneratorApp.decode.opticsSubjectDataBaseDropdown(app, 'valueToSlider', app.opticsParams.subjectDataset);
    
    % The subject ID
    CSFGeneratorApp.decode.opticsSubjectRankSpinner(app, 'valueToSlider', app.opticsParams.subjectRank);
    
    % The pupil size
    CSFGeneratorApp.decode.opticsPupilSizeSpinner(app, 'valueToSlider', app.opticsParams.pupilDiameterMM);
    
    % The wavefront spatial samples
    CSFGeneratorApp.decode.opticsWavefrontSpatialSamplesSpinner(app, 'valueToSlider', app.opticsParams.wavefrontSpatialSamples);
    
    % The zero-center PSF option
    CSFGeneratorApp.decode.opticsZeroCenterPSFCheckBox(app, 'valueToSlider', app.opticsParams.zeroCenterPSF);
    
    % The flipPSFUpsideDown option
    CSFGeneratorApp.decode.opticsFlipPSFUpsideDownCheckBox(app, 'valueToSlider', app.opticsParams.flipPSFUpsideDown);
    
    % The keep optics in sycn with cone mosaic
    CSFGeneratorApp.decode.opticsKeepConeMosaicActivationInSyncCheckBox(app, 'valueToSlider', app.viewModes.opticsKeepConeMosaicActivationInSync);
end

function initializeConeMosaicGUIComponents(app)
    % The view mode - conetypes, cones + retinal image, activation, modulation, or redisual
    CSFGeneratorApp.decode.coneMosaicViewModeKnob(app, 'valueToSlider', app.viewModes.coneMosaic);
    
    % The visualization domain
    CSFGeneratorApp.decode.coneMosaicVisualizationDomainSwitch(app, 'valueToSlider', app.viewModes.coneMosaicVisualizationDomain);

    % The activation type - noise-free or noisy response instances
    CSFGeneratorApp.decode.coneMosaicActivationTypeSwitch(app, 'valueToSlider', app.viewModes.coneMosaicActivationType);
    
    % The activation signal - excitations or photocurrents
    CSFGeneratorApp.decode.coneMosaicActivationSignalSwitch(app, 'valueToSlider', app.viewModes.coneMosaicActivationSignal);

    % The activation dimensionality - 2D space or space-time
    CSFGeneratorApp.decode.coneMosaicActivationDimensionalitySwitch(app, 'valueToSlider', app.viewModes.coneMosaicActivationDimensionality);

    % The ecc-varying options
    CSFGeneratorApp.decode.coneMosaicEccVaryingMacularPigmentDensityCheckBox(app, 'valueToSlider', app.coneMosaicParams.eccVaryingMacularPigmentDensity);
    CSFGeneratorApp.decode.coneMosaicEccVaryingConeApertureCheckBox(app, 'valueToSlider', app.coneMosaicParams.eccVaryingConeAperture);
    CSFGeneratorApp.decode.coneMosaicEccVaryingConeApertureBlurCheckBox(app, 'valueToSlider', app.coneMosaicParams.eccVaryingConeApertureBlur);
    CSFGeneratorApp.decode.coneMosaicEccVaryingOuterSegmentLengthCheckBox(app, 'valueToSlider', app.coneMosaicParams.eccVaryingOuterSegmentLength);
    CSFGeneratorApp.decode.coneMosaicEccVaryingMacularPigmentDynamicCheckBox(app, 'valueToSlider', app.coneMosaicParams.eccVaryingMacularPigmentDynamic);

    % The integration time
    CSFGeneratorApp.decode.coneMosaicIntegrationTimeSpinner(app, 'valueToSlider', app.coneMosaicParams.integrationTimeSeconds);
    
    % The LMS cone ratios
    CSFGeneratorApp.decode.coneMosaicLconeRatioSpinner(app, 'valueToSlider', app.coneMosaicParams.lConeRatio);
    CSFGeneratorApp.decode.coneMosaicMconeRatioSpinner(app, 'valueToSlider', app.coneMosaicParams.mConeRatio);
    CSFGeneratorApp.decode.coneMosaicSconeRatioSpinner(app, 'valueToSlider', app.coneMosaicParams.sConeRatio);
    
    % The tritanopic radius
    CSFGeneratorApp.decode.coneMosaicTritanopicRadiusSpinner(app, 'valueToSlider', app.coneMosaicParams.tritanopicRadiusDegs);
end


function initializeROIGUIComponents(app)
    
    initializeFieldOfView(app);
    initializePolarEccentricity(app);
    if (strcmp(app.roiParams.radialEccentricityScaling, 'log'))
        CSFGeneratorApp.initialize.radialEccentricityWithLogScaling(app);
    else
        CSFGeneratorApp.initialize.radialEccentricityWithLinearScaling(app);
    end
    
    CSFGeneratorApp.decode.roiEyeSwitch(app, 'valueToSlider', app.roiParams.whichEye);
    % Since we are starting with Polans optics, disable the eye switch
    app.roiEyeSwitch.Enable = 'off';
    
    CSFGeneratorApp.decode.roiRadialEccentricitySlider(app, 'valueToSlider', app.roiParams.radialEccentricityDegs);
    CSFGeneratorApp.decode.roiRadialEccentricityScalingSwitch(app, 'valueToSlider', app.roiParams.radialEccentricityScaling);
    CSFGeneratorApp.decode.roiPolarEccentricitySlider(app, 'valueToSlider', app.roiParams.polarEccentricityDegs);
    CSFGeneratorApp.decode.roiFieldOfViewSlider(app, 'valueToSlider', app.roiParams.fieldOfViewDegs);
    
    app.roiMagnificationFactorEditField.Value = app.roiParams.magnificationFactorMicronsPerDeg;

    function initializeFieldOfView(app)
        app.roiFieldOfViewSlider.Limits = [0.0 10];
        app.roiFieldOfViewSlider.MajorTicks = 0:1:10;
        app.roiFieldOfViewSlider.MajorTickLabels = {'0.1', '1', '2', '3', '4', '5', '6','7', '8', '9', '10'};
        app.roiFieldOfViewSlider.MinorTicks = [];
    end

    function initializePolarEccentricity(app)
        app.roiPolarEccentricityKnob.Limits = [0 360];
        app.roiPolarEccentricityKnob.MajorTicks = 0:45:360;
        app.roiPolarEccentricityKnob.MajorTickLabels = {'360', '315', '270', '225', '180', '135',  '90',  '45', '0'};
        app.roiPolarEccentricityKnob.MinorTicks = 0:15:360;

        CSFGeneratorApp.decode.roiPolarEccentricitySlider(app, 'valueToSlider', app.roiParams.polarEccentricityDegs);
    end
end

function initializeComputationalObserverGUIComponents(app)

    % The psychometric function estimation components
    CSFGeneratorApp.decode.psychometricFunctionClassifierTypeDropDown(app, 'valueToSlider', app.psychometricFunctionParams.classifierType);
    CSFGeneratorApp.decode.psychometricFunctionClassifierTrainingTrialsSpinner(app, 'valueToSlider', app.psychometricFunctionParams.trainingTrials);
    CSFGeneratorApp.decode.psychometricFunctionClassifierTestTrialsSpinner(app, 'valueToSlider', app.psychometricFunctionParams.testTrials);
    CSFGeneratorApp.decode.psychometricFunctionContrastLevelsSpinner(app, 'valueToSlider', app.psychometricFunctionParams.contrastLevels);
    CSFGeneratorApp.decode.psychometricFunctionLog10ContrastMinSpinner(app, 'valueToSlider', app.psychometricFunctionParams.log10ContrastMin);
    CSFGeneratorApp.decode.psychometricFunctionLog10ContrastMaxSpinner(app, 'valueToSlider', app.psychometricFunctionParams.log10ContrastMax);
    CSFGeneratorApp.decode.psychometricFunctionLog10ContrastDeltaSpinner(app, 'valueToSlider', app.psychometricFunctionParams.log10ContrastDelta);
    CSFGeneratorApp.decode.psychometricFunctionSlopeMinSpinner(app, 'valueToSlider', app.psychometricFunctionParams.slopeMin);
    CSFGeneratorApp.decode.psychometricFunctionSlopeMaxSpinner(app, 'valueToSlider', app.psychometricFunctionParams.slopeMax);
    CSFGeneratorApp.decode.psychometricFunctionSlopeDeltaSpinner(app, 'valueToSlider', app.psychometricFunctionParams.slopeDelta);
    CSFGeneratorApp.decode.psychometricFunctionEstimationMethodSwitch(app, 'valueToSlider', app.psychometricFunctionParams.estimationMethod);
    
    % The CSF components
    CSFGeneratorApp.decode.csfSpatialFrequencyMinSpinner(app, 'valueToSlider', app.csfParams.spatialFrequencyMin);
    CSFGeneratorApp.decode.csfSpatialFrequencyMaxSpinner(app, 'valueToSlider', app.csfParams.spatialFrequencyMax);
    CSFGeneratorApp.decode.csfSpatialFrequencySamplesSpinner(app, 'valueToSlider', app.csfParams.spatialFrequencySamples);
    CSFGeneratorApp.decode.csfConstantParameterSwitch(app, 'valueToSlider', app.csfParams.constantParameter);
    CSFGeneratorApp.decode.csfNumberOfConstantCyclesSpinner(app, 'valueToSlider', app.csfParams.numberOfConstantCycles);
    CSFGeneratorApp.decode.csfSourceSignalKnob(app, 'valueToSlider', app.csfParams.sourceSignal);
    
end



