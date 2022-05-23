function t_computingWithCustomPSFs()
% Illustrate generate optics using custom PSFs
%
% Syntax:
%   t_computingWithCustomPSFs
%
% Description:
%    Script that demonstrates how to compute cone excitations using optics
%    derived from a set of custom PSFs measured/computed at a set of wavelengths
%    Also illustrates how to compute a text scene realized on a custom display.
% 
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    12/01/21  npc  Wrote it.

    
    % Here we are synthesizing arbitary PSFs.
    % Alternatively we could be importing PSFs measured at a set of wavelengths.
    [thePSFensemble, opticsParams] = synthesizePSFs();
    
    % Generate optics from the synthesized PSFs
    theOI = oiFromPSF(thePSFensemble, opticsParams.wavelengthSupport, ...
        opticsParams.spatialSupportArcMin, opticsParams.pupilDiameterMM, opticsParams.umPerDegree);
    
    visualizeThePSFs = true;
    if (visualizeThePSFs)
        sampledWavelengths = 375:25:700;
        visualizePSFsamples(theOI, sampledWavelengths);
    end
    
    %% Generate a custom experimental display
    % Optional settable key-value pairs are:
    % 'dotsPerInch'                                             - scalar
    % 'viewingDistanceMeters'                                   - scalar
    % 'gammaTable'                                              - [mValues x 3] matrix of LUTs
    % 'wavelengthSupportNanoMeters'                             - [nWaves x 1] matrix of wavelengths
    % 'ambientSPDWattsPerSteradianM2NanoMeter'                  - [nWaves x 1] matrix of the ambient SPD
    % 'spectralPowerDistributionWattsPerSteradianM2NanoMeter'   - [nWaves x 3] matrix of the RGB guns SPDs   
    presentationDisplay = generateCustomDisplay(...
        'dotsPerInch', 220, ...
        'viewingDistanceMeters', 1.00, ...
        'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
        'plotCharacteristics', true);
    
    %% Select stimulus chromaticity specification
    chromaSpecificationType = 'RGBsettings';  % choose between {'RGBsettings', 'chromaLumaLMScontrasts'}
    
    % gamma uncorrected values
    switch (chromaSpecificationType)
        case 'RGBsettings'
        % Specify both background and stimulus in terms of RGB settings
        % values on presentation display
        chromaSpecification = struct(...
                'type', chromaSpecificationType, ...
                'backgroundRGB', [0.5 0.5 0.5], ...
                'foregroundRGB',  [1 1 1]);
                    
        case 'chromaLumaLMScontrasts'
        % Specify background in terms of cie-31 (x,y) and luminance (cd/m2)
        % and stimulus in terms of nominal LMS cone contrasts
        chromaSpecification = struct(...
               'type', chromaSpecificationType, ...
               'backgroundChromaLuma', [0.31 0.32 40], ...
               'foregroundLMSConeContrasts', [-0.5 -0.5 0.0]); 
    end
    
    % Stimulus params
    theString = 'Hello ISETBio world ! ';
    textSceneParams = struct(...
        'textString', theString, ...                   % Text to display
        'textRotation', 0, ...                         % Rotation (0,90,180,270 only)
        'rowsNum', 60, ...                             % Pixels along the vertical (y) dimension
        'colsNum', 400, ...                            % Pixels along the horizontal (x) dimension
        'targetRow', 20, ...                           % Stimulus Y-pixel offset 
        'targetCol', 20, ...                           % Stimulus X-pixel offset 
        'chromaSpecification', chromaSpecification ... % Background and stimulus chromaticity
    );
            
    % Generate the scene
    visualizeScene = true;
    theScene = rotatedTextSceneRealizedOnDisplay(presentationDisplay, textSceneParams, visualizeScene);
    
    %% Compute the optical image
    theOI = oiCompute(theOI, theScene);
    
    % Visualize the optical image
    visualizeOpticalImage(theOI, 'crossHairsAtOrigin', true, 'displayRadianceMaps', false);
    
    %% Generate the cone mosaic 
    % Specify a size matching the scene and a retinal
    % magnification factor, matching that of the optics
    optics = oiGet(theOI, 'optics');
    focalLengthMicrons = opticsGet(optics, 'focal length')*1e6;
    micronsPerDegree = focalLengthMicrons*tand(1);
    
    % Generate a cone mosaic whose size matches the size of the stimulus
    fovDegreesWidth = sceneGet(theScene, 'wAngular');
    fovDegreesHeight = sceneGet(theScene, 'hAngular');
    theConeMosaic = cMosaic(...
        'sizeDegs', [fovDegreesWidth fovDegreesHeight], ...
        'micronsPerDegree', micronsPerDegree, ...
        'integrationTime', 2/1000);
    
    % Generate 1 eye movement sequence lasting for 100 msec
    eyeMovementDurationSeconds = 200/1000;
    theConeMosaic.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'nTrials', 1, ...
        'randomSeed', 10);
    
    % Visualize cone mosaic and eye movement path
    theConeMosaic.visualize('displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:100));
    
    
    %% Compute the cone mosaic response time series
    [coneExcitationsNoiseFree, coneExcitations] = ...
        theConeMosaic.compute(theOI); %, 'withFixationalEyeMovements', true);
    
    % Visualize the scene, the optical image and the cone mosaic activation together
    visualizeEverything(theScene, theOI, theConeMosaic, coneExcitations, ...
        fovDegreesWidth, fovDegreesHeight);
end




function [thePSFensemble, opticsParams] = synthesizePSFs()

    % Optics params 
    opticsParams = struct(...
        'wavelengthSupport', 400:5:700, ...     % list of wavelengths at which we have measured/synthesized PSF data
        'spatialSupportArcMin', -10:0.1:10, ... % list of spatial samples at which we have measured/synthesized PSF data
        'pupilDiameterMM', 3.0, ...             % pupil diameter at which PSFs were measured
        'umPerDegree', 292 ...                  % microns-per-degree retinal magnification
        );
    
    % Synthesized PSFs have a 90 deg at the in-focus wavenength
    rotation0 = 90;
    
    % Synthesized PSFs are in-focus at 560 nm
    inFocusWavelength = 560;
    
    % In-focus synthesized PSF has a sigma of 0.5 arc min
    inFocusSigmaArcMin = 0.15;
   
    % Synthesized PSFs shift with wavelength with a rate of 0.007 arc min/nm 
    lateralShiftArcMinPerNM = 0.01;
    
    % Synthesized PSFs change their orientation with this rate / nm
    rotationChangePerNM = 45/150;
    
    % Synthesized PSFs become progressively defocused from the in-focus
    % wavelength with a rate of 0.004 arc min/nm
    sigmaIncreaseArcMinPerNM = 0.004;
    
    % Spatial grid
    [x,y] = meshgrid(opticsParams.spatialSupportArcMin,opticsParams.spatialSupportArcMin);
    
    % Preallocate memory
    thePSFensemble = zeros(numel(opticsParams.spatialSupportArcMin), numel(opticsParams.spatialSupportArcMin), numel(opticsParams.wavelengthSupport));
    
    for wIndex = 1:numel(opticsParams.wavelengthSupport)
        % simuate lateral shift with wavelength
        deltaWavelength = opticsParams.wavelengthSupport(wIndex)-inFocusWavelength;
        xo = 0 + deltaWavelength*lateralShiftArcMinPerNM;
        yo = 0;
        
        % simulate defocus with wavelength
        sigmaArcMin = inFocusSigmaArcMin + abs(deltaWavelength)*sigmaIncreaseArcMinPerNM;
        
        % rotate axes
        rotation = rotation0 - deltaWavelength*rotationChangePerNM;
        xx = (x-xo) * cosd(rotation) + (y-yo) * sind(rotation);
        yy =-(x-xo) * sind(rotation) + (y-yo) * cosd(rotation);
        
        % Generate model PSF
        thePSFensemble(:,:,wIndex) = ...
                 exp(-0.5*(xx/sigmaArcMin/2).^2) .* ...
                 exp(-0.5*(yy/sigmaArcMin).^2);
    end
    
end

function visualizeEverything(theScene, theOI, theConeMosaic, coneExcitations, ...
    fovDegreesWidth, fovDegreesHeight)

    % Visualize scene and optical image
    hFig = figure(100); clf;
    set(hFig, 'Position', [10 10 2400 1300], 'Color', [ 1 1 1]);
    colsNum = 1;
    rowsNum = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rowsNum, ...
               'colsNum', colsNum, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.03, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.035, ...
               'topMargin',      0.02);
           
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualizeScene(theScene, ...
            'axesHandle', ax, ...
            'spatialSupportInDegs', true, ...
            'crossHairsAtOrigin', true, ...
            'displayRadianceMaps', false);
    set(ax, 'XLim', 0.5*fovDegreesWidth*[-1 1], 'YLim', 0.5*fovDegreesHeight*[-1 1]);
    drawnow;
    
    ax = subplot('Position', subplotPosVectors(2,1).v);
    visualizeOpticalImage(theOI, 'axesHandle', ax, 'crossHairsAtOrigin', true, 'displayRadianceMaps', false);
    set(ax, 'XLim', 0.5*fovDegreesWidth*[-1 1], 'YLim', 0.5*fovDegreesHeight*[-1 1]);
    drawnow;
    
    ax = subplot('Position', subplotPosVectors(3,1).v);

    for iTrial = 1:size(coneExcitations,1)
        for tBin = 1:size(coneExcitations,2)
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'activation', squeeze(coneExcitations(iTrial, tBin,:)), ...
            'activationRange', [min(coneExcitations(:)) max(coneExcitations(:))], ...
            'domainVisualizationLimits', [-0.5*fovDegreesWidth 0.5*fovDegreesWidth -0.5*fovDegreesHeight 0.5*fovDegreesHeight], ...
            'domainVisualizationTicks', struct('x', -2.5:0.5:2.5, 'y', -2.5:0.5:2.5), ...
            'backgroundColor', [0 0 0], ...
            'verticalActivationColorBarInside', true, ...
            'crossHairsOnMosaicCenter', true, ...
            'crossHairsColor', [0 0 0], ...
            'plotTitle', sprintf('the cone mosaic activation (time: %2.0f msec)', ...
                (tBin-1) * theConeMosaic.integrationTime*1000));
        drawnow;
        end
    end
end



function visualizePSFsamples(theOI, sampledWavelengths)
    hFig = figure(); clf;
    set(hFig,  'Color', [1 1 1], 'Position', [10 10 1800 600]);
    colsNum = 7;
    rowsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rowsNum, ...
               'colsNum', colsNum, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.05, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);
        
     psfRangeArcMin = 8;
     wavelengthSupport = oiGet(theOI, 'wave');
     for k = 1:numel(sampledWavelengths)
        [~,wIndex] = min(abs(wavelengthSupport-sampledWavelengths(k)));
        targetWavelength = wavelengthSupport(wIndex);
        
        row = floor((k-1)/colsNum)+1;
        col = mod(k-1,colsNum)+1;
        ax = subplot('Position', subplotPosVectors(row,col).v);
        visualizePSF(theOI, targetWavelength, psfRangeArcMin, ...
            'contourLevels', 0.1:0.1:0.9, ...
            'axesHandle', ax, ...
            'figureTitle', sprintf('%2.0f nm', targetWavelength), ...
            'fontSize', 14);
     end
end