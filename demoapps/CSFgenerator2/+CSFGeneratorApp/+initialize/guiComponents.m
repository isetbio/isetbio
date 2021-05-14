function guiComponents(app)

    % Initialize colors for different aspects of the app
    initializeColors(app)
    
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
    app.statusMessages('performance assessment') = s;
    app.statusMessages('optics') = s;
    app.statusMessages('cone mosaic') = s;
    app.statusMessages('fixational eye movements') = s;
    app.statusMessages('mRGC mosaic') = s;
    app.statusMessages('pRGC mosaic') = s;
    
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
    
    % Stimulus spatial envelope, orientation, and spatial frequency
    CSFGeneratorApp.decode.stimulusSpatialEnvelopeDropDown(app, 'valueToSlider', app.stimParams.spatialEnvelope);
    CSFGeneratorApp.decode.stimulusSpatialFrequencySpinner(app, 'valueToSlider', app.stimParams.spatialFrequencyCPD);
    CSFGeneratorApp.decode.stimulusOrientationSpinner(app, 'valueToSlider', app.stimParams.orientationDegs);
end

function initializeOpticsGUIComponents(app)
    % The visualized wavelength slider
    app.opticsVisualizedWavelengthSlider.Limits = [400 750];
    app.opticsVisualizedWavelengthSlider.MajorTicks = 400:50:750;
    app.opticsVisualizedWavelengthSlider.MinorTicks = 400:10:750;
    CSFGeneratorApp.decode.opticsWavelengthSlider(app, 'valueToSlider', app.opticsParams.visualizedWavelength);
    
    % The subject dataset
    CSFGeneratorApp.decode.opticsSubjectDataBaseDropdown(app, 'valueToSlider', app.opticsParams.subjectDataset);
    
    % The subject ID
    CSFGeneratorApp.decode.opticsSubjectIDSpinner(app, 'valueToSlider', app.opticsParams.subjectID);
    
    % The pupil size
    CSFGeneratorApp.decode.opticsPupilSizeSpinner(app, 'valueToSlider', app.opticsParams.pupilDiameterMM);
    
    % The central refraction subtraction
    CSFGeneratorApp.decode.opticsSubtractCentralRefractionCheckBox(app, 'valueToSlider', app.opticsParams.subtractCentralRefraction);
    
    % The central refraction subtraction
    CSFGeneratorApp.decode.opticsKeepConeMosaicActivationInSyncCheckBox(app, 'valueToSlider', app.viewModes.opticsKeepConeMosaicActivationInSync);
end

function initializeConeMosaicGUIComponents(app)
    % The view mode
    CSFGeneratorApp.decode.coneMosaicViewModeKnob(app, 'valueToSlider', app.viewModes.coneMosaic);
    
    % The activation type
    CSFGeneratorApp.decode.coneMosaicActivationTypeSwitch(app, 'valueToSlider', app.viewModes.coneMosaicActivationType);
    
    % The visualization domain
    CSFGeneratorApp.decode.coneMosaicVisualizationDomainSwitch(app, 'valueToSlider', app.viewModes.coneMosaicVisualizationDomain);
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


