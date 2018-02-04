function t_fixationalEyeMovementsToPhotocurrents

    % Generate oiSequence
    fovDegs = 0.1;
    theOIsequence = generateImpulseOIsequence(fovDegs);

    % Instantiate a cone mosaic
    integrationTime = 5/1000;
    
    cm = coneMosaicGenerate(fovDegs,integrationTime);
    
    % Instantiate a fixational eye movement object
    fixEMobj = fixationalEM();
    
    % Generate emPaths for this mosaic
    % (we need to know the mosaic we are generating emPaths for
    % so as to obtain the integrationTime, the patternSampleSize 
    % and the micronsPerDegree conversion factor)
    nTrials = 2;
    eyeMovementsPerTrial = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
    fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, 'nTrials', nTrials, 'rSeed', 1);
    
    visualizedTrial = 1;
    eyeMovementsData = struct(...
         'show', true, ...
         'timeAxisMillisecs', fixEMobj.timeAxis*1000, ...
         'posMicrons', squeeze(fixEMobj.emPosMicrons(visualizedTrial,:,:)) ...
    );
    figure();
 	theOIsequence.visualize('montage', ...
         'showIlluminanceMap', false, ...
         'eyeMovementsData', eyeMovementsData);
 

    % Visualize one trial of eye movements on top of the mosaic.
    % Here we use the 'emPosArcMin' data, which contain the contain the eye movement paths in units of patternSampleSize
    
    figure(1); 
    subplot(1,3,1);
    hold on
    plot(fixEMobj.timeAxis*1000, squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,1)), 'r-', 'MarkerSize', 5);
    subplot(1,3,2);
    hold on
    plot(fixEMobj.timeAxis*1000, squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,2)), 'b-', 'MarkerSize', 5);
    subplot(1,3,3);
    hold on
    velocityMeasurementIntervalSeconds = 20/1000;
    velocity = fixEMobj.computeVelocity(fixEMobj.timeAxis, squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,:)), velocityMeasurementIntervalSeconds);
    plot(fixEMobj.timeAxis*1000, velocity, 'k-', 'MarkerSize', 5);
    
    % Visualize one trial of eye movements on top of the mosaic.
    % Here we use the 'emPosMicrons' data, which contain the contain the eye movement paths in units of patternSampleSize
    
    
    figure(10);
    subplot(1,2,1)
    cm.visualizeGrid(...
        'axesHandle', gca, ...
        'overlayEMpathMicrons', squeeze(fixEMobj.emPosMicrons(visualizedTrial,:,:)), ...
        'overlayNullSensors', false, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', false,...
        'generateNewFigure', false);
    
    
    
    % Compute the mosaic response. 
    % Here we use the 'emPos' data, which contain the eye movement paths in units of patternSampleSize
    fprintf('Compute mosaic cone responses (isomerizations and photocurrents) to the background sequence, taking eye movements into account.\n');
    [isomerizationsAdaptingStim, photocurrentsAdaptingStim, osImpulseResponseFunctions, osMeanCurrents] = ...
    cm.computeForOISequence(theOIsequence, ...
        'emPaths', fixEMobj.emPos, ...
        'interpFilters', [], ...
        'meanCur', [], ...
        'currentFlag', true);
    
    if (1==2)
        mosaicExtentMicrons = 0.5*max(cm.fov(:))*cm.micronsPerDegree;
        emPathsRangeMicrons = max(abs(fixEMobj.emPosMicrons(:)));
        spatialExtent = max([mosaicExtentMicrons emPathsRangeMicrons ]) * 1.05 * [-1 1];
        ticksMeters = (-100:10:100)*1e-6;
        for k = 1:numel(ticksMeters)
            if (mod(k,2) == 1)
                tickLabels{k} = sprintf('%2.0f', ticksMeters(k)*1e6);
            else
                tickLabels{k} = ' ';
            end
        end
        set(gca, 'XTick', ticksMeters, 'YTick', ticksMeters, 'XTickLabel', tickLabels, 'YTickLabel', tickLabels);
        set(gca, 'XLim', spatialExtent*1e-6, 'YLim', spatialExtent*1e-6);
    end
    
end

function cm = coneMosaicGenerate(fovDegs, integrationTime)
    resamplingFactor = 11;
    eccBasedConeDensity = false;
    customLamda = [];
    customInnerSegmentDiameter = [];
    
    % Instantiate a hex mosaic
    cm = coneMosaicHex(resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'eccBasedConeDensity', eccBasedConeDensity, ... 
        'customLambda', customLamda, ...
        'customInnerSegmentDiameter', customInnerSegmentDiameter ...
    );

    % Set the mosaic's integration time
    cm.integrationTime = integrationTime;
end

function theOIsequence = generateImpulseOIsequence(fovDegs)
    sparams.fov = fovDegs; sparams.luminance = 100;
    stimWeights = zeros(1,50); 
    stimWeights(4) = 1;
    stimRefreshInterval = 20/1000;
    
    theOIsequence = oisCreate('impulse','add', stimWeights, ...
        'sampleTimes', stimRefreshInterval*((1:length(stimWeights))-1), ...
        'sceneParameters',sparams);
    %figure();
    %theOIsequence.visualize('format', 'montage', 'showIlluminanceMap', true);   
end
