function t_osTimeStep
%%t_osTimeStep  Show how to change timebase of os calculation.
%
% Description:
%    Demonstrate simulations using three different timebases, one for
%    stimuli (based on stimulus refresh rate), one for absorptions and eye
%    movements (based on coneMosaic.integrationTime), and a third one for
%    outer segment current computations (based on os.timeStep)
%
%    Also demonstrates usage of the computeForOISequence() method of
%    @coneMosaic, which computes absorptions and photocurrents for a
%    sequence of sequentially presented optical images with eye movements.
%
%    This is a rather complex tutorial that can demonstrate lots of things.
%    It has many subfunctions which do lots of the work, and which have no
%    comments.
%
%    [DHB NOTE: This needs more extensive comments before it can be
%    used in an introductory way. For NPC to take on.]

% NPC, ISETBIO TEAM, 2016
%
% 08/07/17  dhb  Commenting pass.  Decided this was too hard and moved to advanced.


%% Init
ieInit;

%% Examine the effects of varying the integrationTime
conditionSet = 1;

% Examine effects of varying the response time interval (os.timeStep)
%conditionSet = 2;

% Examine and contrast the magnitudes of photon noise vs OS noise
%conditionSet = 3;

% Custom condition
%conditionSet = 4;

condData = makeConditionSet(conditionSet);

% Run all the conditions
for stimulusConditionIndex = 1:numel(condData)
    % Get the condition data
    c = condData{stimulusConditionIndex};   
    
    % Run the simulation for this condition
    [theConeMosaic, theOIsequence, ...
        isomerizationRateSequence, photoCurrentSequence, ...
        oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis] = runSimulation(c.mosaicSize, c.meanLuminance, c.modulation, c.modulationRegion, c.stimulusSamplingInterval, c.integrationTime, c.osTimeStep, c.photonNoise, c.osNoise);
    
    % Plot the results
    plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis, stimulusConditionIndex, c);
end

end

function [theConeMosaic, theOIsequence, ...
    isomerizationRateSequence, photoCurrentSequence, ...
    oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis] = runSimulation(mosaicSize, meanLuminance,  modulationGain, modulationRegion, stimulusSamplingInterval, integrationTime, osTimeStep, photonNoise, osNoise)

    % Define the time axis for the simulation
    minTime = -0.84;
    maxTime = 0.72;
    oiTimeAxis = minTime:stimulusSamplingInterval:maxTime;
    
    % Compute the stimulus modulation function
    stimulusRampTau = 0.18;
    modulationFunction = modulationGain * exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);

    % Generate a uniform field scene with desired mean luminance
    if (isnan(mosaicSize))
        FOV = 0.2;
    else
        FOV = max(mosaicSize);
    end
    theScene = uniformFieldSceneCreate(FOV, meanLuminance);

    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);

    % Generate the sequence of optical images
    theOIsequence = oiSequenceGenerate(theScene, theOI, oiTimeAxis, modulationFunction, modulationRegion);

    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep);

    % Generate eye movement sequence for all oi's
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    theConeMosaic.emGenSequence(eyeMovementsNum);
    
    [absorptionsCountSequence, photoCurrentSequence] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'currentFlag', true ...
            );
    absorptionsTimeAxis = theConeMosaic.timeAxis + theOIsequence.timeAxis(1);  
    photoCurrentTimeAxis = absorptionsTimeAxis;
    
    %theOIsequence.visualizeWithEyeMovementSequence(absorptionsTimeAxis);
    
    % Compute photon rate from photon count
    isomerizationRateSequence = absorptionsCountSequence / theConeMosaic.integrationTime;
end


% ------- Helper functions --------

function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicSize)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicSize);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;
end


function theOIsequence = oiSequenceGenerate(theScene, theOI, oiTimeAxis, modulationFunction, modulationType)
    % Compute the background and modulated optical images
    oiBackground = oiCompute(theOI, theScene);
    oiModulated  = oiBackground;
    
    if strcmp(modulationType, 'FULL')
        theOIsequence = oiSequence(oiBackground, oiModulated, oiTimeAxis, modulationFunction, 'composition', 'add');
    else
        pos = oiGet(oiBackground, 'spatial support', 'microns');
        modulationRegion.radiusInMicrons = 0.75*max(pos(:));
        theOIsequence = oiSequence(oiBackground, oiModulated, oiTimeAxis, modulationFunction, 'modulationRegion', modulationRegion);
    end
end


function theOI = oiGenerate(noOptics)
    % Generate optics
    if (noOptics)
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('human');
    end
end


function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniform equal photon', 128);
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end


function plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, oiTimeAxis, absorptionsTimeAxis, responseTimeAxis, figNo, condData)

    % Plot the sequence of OIs with the eye movements 
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10+figNo*50 10+figNo*100 1920 760], 'Color', [1 1 1]);
    set(hFig, 'Name', sprintf('Scene Mean Luminance: %2.1f cd/m2,     Modulation: %2.2f,     Stimulus Sampling: %2.1f ms,     Integration Time: %2.1f ms,   osTimeStep: %2.1f ms,      PhotonNoise: %s,      osNoise: %s', condData.meanLuminance, condData.modulation, condData.stimulusSamplingInterval*1000, condData.integrationTime*1000, condData.osTimeStep*1000, condData.photonNoise, condData.osNoise));

    tabGroup = uitabgroup('Parent', hFig);
    
    eyeMovementsTab = uitab(tabGroup, 'Title', '2D eye movements and OI sequence', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 1]);
    timeSeriesTab = uitab(tabGroup, 'Title', 'Time series: optical image photon rate, eye movements, aborptions, & photocurrents', 'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 1]);
    
    set(tabGroup, 'SelectedTab',eyeMovementsTab);
    axes('parent',eyeMovementsTab);

    
    plotRows = round(0.75*sqrt(theOIsequence.length)); 
    plotCols = ceil(theOIsequence.length/plotRows);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', plotRows, ...
           'colsNum', plotCols, ...
           'heightMargin',   0.03, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02);
       
    maxRGB = 0;
    for oiIndex = 1:theOIsequence.length
        tmp = xyz2rgb(oiGet(theOIsequence.frameAtIndex(oiIndex), 'xyz'));
        if (maxRGB < max(tmp(:)))
            maxRGB = max(tmp(:));
        end
         oiImage{oiIndex} = tmp;
    end
    
    % Retrieve the eye movement sequence
    eyeMovementSequence = theConeMosaic.emPositions;
    
    for oiIndex = 1:theOIsequence.length
        pos = oiGet(theOIsequence.frameAtIndex(oiIndex), 'spatial support', 'microns');
        oiXaxis = pos(1,:,1); oiYaxis = pos(:,1,2);
        r = floor((oiIndex-1)/plotCols)+1;
        c = mod((oiIndex-1), plotCols)+1;
        subplot('Position', subplotPosVectors(r,c).v);
        
        % Plot the OI
        imagesc(oiXaxis, oiYaxis, oiImage{oiIndex}/maxRGB);
        % Overlay the eye movements up to this point
        hold on;
        idx = find(absorptionsTimeAxis < oiTimeAxis(oiIndex));
        plot(eyeMovementSequence(idx,1)*theConeMosaic.pigment.width*1e6, eyeMovementSequence(idx,2)*theConeMosaic.pigment.width*1e6, 'ks-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.3 0.3 0.3]);
        if (oiIndex < theOIsequence.length)
            idx = find((absorptionsTimeAxis>=oiTimeAxis(oiIndex)) & (absorptionsTimeAxis<oiTimeAxis(oiIndex+1)));
        else
            idx = find((absorptionsTimeAxis>=oiTimeAxis(oiIndex)));
        end
        % Emphasize in red, the eye movements for the current framer
        plot(eyeMovementSequence(idx,1)*theConeMosaic.pigment.width*1e6, eyeMovementSequence(idx,2)*theConeMosaic.pigment.width*1e6, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', [1.0 0.5 0.5]);
        
        % overlay the cone mosaic
        if (oiIndex == 1)
            plot(theConeMosaic.coneLocs(:,1)*1e6, theConeMosaic.coneLocs(:,2)*1e6, 'k.');
        end
        
        axis 'image'; axis 'xy';
        set(gca, 'CLim', [0 1], 'FontSize', 14);
        if ~((r == plotRows) && (c == 1))
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('t: %2.1f msec', 1000*oiTimeAxis(oiIndex)), 'FontSize', 14);
    end
    
    % Plot time-series responses
    set(tabGroup, 'SelectedTab',timeSeriesTab);
    axes('parent',timeSeriesTab);
    
    %% Plot the photon rate at the center of the optical image
    subplot('Position', [0.03 0.07 0.18 0.89]);
    oiWavelengthAxis = oiGet(theOIsequence.frameAtIndex(1), 'wave');
    referencePositionOpticalImagePhotons = zeros(numel(oiWavelengthAxis), theOIsequence.length);
    for oiIndex = 1:theOIsequence.length
        retinalPhotonsAtCurrentFrame = oiGet(theOIsequence.frameAtIndex(oiIndex), 'photons');
        refRow = round(size(retinalPhotonsAtCurrentFrame,1)/2);
        refCol = round(size(retinalPhotonsAtCurrentFrame,2)/2);
        referencePositionOpticalImagePhotons(:, oiIndex) = squeeze(retinalPhotonsAtCurrentFrame(refRow, refCol, :));
    end
    hP = pcolor(oiTimeAxis, oiWavelengthAxis, referencePositionOpticalImagePhotons);
    set(hP, 'EdgeColor', 'none');
    hold on;
    % Plot the total photons (summed across all wavelengths)
    totalPhotons = sum(referencePositionOpticalImagePhotons,1);
    totalPhotonsNorm = oiWavelengthAxis(1) + (oiWavelengthAxis(end)-oiWavelengthAxis(1))*(totalPhotons-min(totalPhotons))/(max(totalPhotons)-min(totalPhotons));
    stairs(oiTimeAxis, totalPhotonsNorm, 'c-', 'LineWidth', 2.0);
    plotStimulusTimes([min(totalPhotonsNorm) max(totalPhotonsNorm)]);
    
    hold off; box on
    set(gca, 'YLim', [oiWavelengthAxis(1) oiWavelengthAxis(end)], 'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('wavelength (nm)', 'FontSize', 14, 'FontWeight', 'bold');
    title('optical image photon rate (image center)', 'FontSize', 14);
    hC = colorbar('Location', 'NorthOutside', 'parent', timeSeriesTab);
    hC.FontSize =  14;
    hC.Label.String = 'photons/sec';
    axis 'xy'
    colormap(gray(1024));

    %% Plot the eye movement sequence (different colors for different OIs)
    subplot('Position', [0.25 0.07 0.22 0.89]); hold on;
    eyeMovementRange = [-100 100];

    plot(absorptionsTimeAxis, eyeMovementSequence(:,1)*theConeMosaic.pigment.width*1e6, '.', 'MarkerSize', 15, 'Color', 'r');
    hold on;
    plot(absorptionsTimeAxis, eyeMovementSequence(:,2)*theConeMosaic.pigment.height*1e6, '.', 'MarkerSize', 15, 'Color', 'b');
    plotStimulusTimes(eyeMovementRange);
    
    box on
    set(gca, 'YLim', [eyeMovementRange(1) eyeMovementRange(end)], 'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    legend({'eye position (X)', 'eye position (Y)'});
    ylabel('X,Y eye position (microns)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('X,Y eye movements', 'FontSize', 14);
    
    
    %% Plot the LMS isomerizations
    if (theConeMosaic.rows ==1) && (theConeMosaic.cols == 3)
       referenceConeRows = [1 1 1]; referenceConeCols = [1 2 3];
    else
       % Find the (row,col) coords of the center-most L, M and S-cone
       for k = 1:3
            coneIndices = find(theConeMosaic.pattern == k+1); 
            [~, idx] = min(sum((theConeMosaic.coneLocs(coneIndices,:)).^2, 2));
            [referenceConeRows(k), referenceConeCols(k)] = ind2sub(size(theConeMosaic.pattern), coneIndices(idx));
       end
    end
    subplot('Position', [0.50 0.07 0.22 0.89]);
    isomerizationRange = [min(isomerizationRateSequence(:)) 1.05*max(isomerizationRateSequence(:))];
    hold  on
    coneColors = [1 0 0; 0 1 0; 0 0 1];
    for k = 1:3
        plot(absorptionsTimeAxis, squeeze(isomerizationRateSequence(1, referenceConeRows(k),referenceConeCols(k),:)), '.', 'Color', squeeze(coneColors(k,:)), 'MarkerSize', 15, 'LineWidth', 1.5);
    end
    plotStimulusTimes(isomerizationRange);
    
    hold off;
    set(gca, 'YLim', isomerizationRange, 'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    ylabel('isomerization rate (R*/cone/sec)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('L,M,S-cone isomerization rates', 'FontSize', 14);

    %% Plot the photocurrents
    subplot('Position', [0.75 0.07 0.22 0.89]);
    photoCurrentRange = [min(photoCurrentSequence(:)) max(photoCurrentSequence(:))+2];
    hold on;
    for k = 1:3
        plot(responseTimeAxis, squeeze(photoCurrentSequence(1, referenceConeRows(k),referenceConeCols(k),:)), 'k.', 'Color', squeeze(coneColors(k,:)), 'MarkerSize', 15, 'LineWidth', 1.5);
    end
    plotStimulusTimes(photoCurrentRange);
    
    hold off;
    set(gca, 'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'YLim', photoCurrentRange, 'FontSize', 14);
    ylabel('photocurrent (pA)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('@osLinear response', 'FontSize', 14);
    
    % Switch to eye movements tab
    set(tabGroup, 'SelectedTab',eyeMovementsTab);
    
    drawnow
    
    function plotStimulusTimes(signalRange)
        % Plot lines demarkating each OI time duration
        for oiIndex = 1:theOIsequence.length
            plot(oiTimeAxis(oiIndex)*[1 1], signalRange, 'k-');
        end
        % Plot the origin in magenta
        plot([0 0], signalRange, '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    end

end


function condData = makeConditionSet(conditionSet)

    % scene mean luminance
    meanLuminance = 1500;    
    
    % Assemble conditions to run.
    condData = {};

    switch conditionSet

        % Effects of varying the integrationTime
        case 1
            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % 50%  modulation against background
                'modulationRegion', 'CENTER', ...           % modulate the central image (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/50, ...      % 50 Hz stimulus refresh, e.g., 20 msec per optical image
                'osTimeStep', 0.1/1000, ...                 % 0.1 millisecond
                'integrationTime', nan, ...                 % we will vary this one
                'photonNoise', 'random', ...                % select from {'random', 'frozen', 'none'}
                'osNoise', 'none');                         % select from {'random', 'frozen', 'none'}
            
            % Varied params
            c0.integrationTime = 100/1000;                  % 100 ms (longer than the stimulus sampling interval, so < 1 eye movement / oi)
            condData{numel(condData)+1} = c0;
            
            c0.integrationTime = 10/1000;                   % 10 ms (smaller than the stimulus sampling interval, 2 eye movements / oi)
            condData{numel(condData)+1} = c0;
            
            c0.integrationTime = 1/1000;                    % 1 ms (20 eye movements / oi)
            condData{numel(condData)+1} = c0;
            
        % Examine effects of varying the os.timeStep
        case 2
            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % 50%  modulation against background
                'modulationRegion', 'CENTER', ...           % modulate the central image (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/10, ...      % 100 Hz stimulus refresh
                'integrationTime', 20/1000, ...             % 20 milliseconds
                'osTimeStep', nan, ...                      % we'll vary that one
                'photonNoise', 'random', ...                    % select from {'random', 'frozen', 'none'}
                'osNoise', 'none');                         % select from {'random', 'frozen', 'none'}

            % Varied params
            c0.osTimeStep = 1/1000;                    % 1 milliseconds
            condData{numel(condData)+1} = c0;
            
            c0.osTimeStep = 0.5/1000;                  % 0.5 milliseconds
            condData{numel(condData)+1} = c0;
            
            c0.osTimeStep = 0.1/1000;                  % 0.1 milliseconds
            condData{numel(condData)+1} = c0;
            
        
            
        % Examine and contrast the magnitudes of photon noise vs OS noise
        case 3    
         
            stimulusRefreshRateInHz = 25;
            eyeMovementsPerStimulusRefresh = 6;
            
            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % % modulation against background
                'modulationRegion', 'CENTER', ...           % modulate the center only  (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/stimulusRefreshRateInHz, ...      
                'integrationTime', 5/1000, ...             % 5 milliseconds
                'osTimeStep', 0.1/1000, ...  
                'photonNoise', nan, ...
                'osNoise', nan);
            
            % Varied params
            % No noise
            c0.photonNoise = 'none';        % select from {'random', 'frozen', 'none'}
            c0.osNoise = 'none';            % select from {'random', 'frozen', 'none'}
            condData{numel(condData)+1} = c0;
            
            % Photon noise
            c0.photonNoise = 'random';      % select from {'random', 'frozen', 'none'}
            c0.osNoise = 'none';            % select from {'random', 'frozen', 'none'}
            condData{numel(condData)+1} = c0;
            
            % OS noise 
            c0.photonNoise = 'none';        % select from {'random', 'frozen', 'none'}
            c0.osNoise = 'random';          % select from {'random', 'frozen', 'none'}
            condData{numel(condData)+1} = c0;
            
            % Both photon noise and OS noise
            c0.photonNoise = 'random';      % select from {'random', 'frozen', 'none'}
            c0.osNoise = 'random';          % select from {'random', 'frozen', 'none'}
            condData{numel(condData)+1} = c0;
            
        % Custom condition - customize it !
        case 4
            
            % Custom
            stimulusRefreshRateInHz = 66;
            integrationTime = 50/1000;
            
            % Steady params
            c0 = struct(...
                'mosaicSize', 0.5, ...                                           % FOV = 0.5 x 0.5 deg
                'meanLuminance', meanLuminance, ...                              % scene mean luminance
                'modulation', 0.5, ...                                           % modulation against background
                'modulationRegion', 'CENTER', ...                                % modulate the center only  (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/stimulusRefreshRateInHz, ...      
                'integrationTime', integrationTime, ...     
                'osTimeStep', 0.1/1000, ...                                      % 0.1 millisecond
                'photonNoise', 'none', ...                                       % select from {'random', 'frozen', 'none'}
                'osNoise', 'none');                                              % select from {'random', 'frozen', 'none'}
            
           condData{numel(condData)+1} = c0;
    end
end
