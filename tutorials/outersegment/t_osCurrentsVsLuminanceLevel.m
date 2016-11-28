function t_osCurrentsVsLuminanceLevel

    % Define the time axis for the simulation
    stimulusSamplingInterval = 50/1000;
    oiTimeAxis = 0:stimulusSamplingInterval:2.0;
    oiTimeAxis = oiTimeAxis - mean(oiTimeAxis);
    
    % Compute the stimulus modulation function
    stimulusRampTau = 0.18;
    modulationFunction = exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);
    
    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);
    
    % Luminance levels to be examined
    luminancesExamined = [100 500 1000 1500 2000 3000 4000 8000];

    % outer segment time steps to be examined
    osTimeSteps = [0.001 0.01 0.1 0.3 0.5 1 3 5]/1000;

    % Generate a number of oiSequences, one for each luminance examined
    for lumIndex = 1:numel(luminancesExamined)
        % Generate the scene
        FOV = 1.0;
        theScene = uniformFieldSceneCreate(FOV, luminancesExamined(lumIndex));

        % Compute the optical images
        oiBackground = oiCompute(theOI, theScene);
        oiModulated  = oiBackground;
    
        % Generate the sequence of optical images
        pos = oiGet(oiBackground, 'spatial support', 'microns');
        modulationRegion.radiusInMicrons = 0.5*max(pos(:));
        theOIsequenceArray{lumIndex} = oiSequence(oiBackground, oiModulated, oiTimeAxis, modulationFunction, 'modulationRegion', modulationRegion);
        %theOIsequence.visualize('format', 'montage');
    end % lumIndex
    
    
    % Generate a cone mosaic with 1 L-, 1 M-, and 1 S-cone only
    mosaicSize =  nan;                    
    integrationTime = 10/1000;            
    photonNoise = false;                     
    osTimeStep = 1/1000;                    
    osNoise = false;                        
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep);
    
    % Generate eye movement path
    eyeMovementsNum = theOIsequenceArray{1}.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    instancesNum = 1;
    emPaths = zeros(instancesNum, eyeMovementsNum,2);
    for instanceIndex = 1:instancesNum
        emPaths(instanceIndex, :, :) = theConeMosaic.emGenSequence(eyeMovementsNum)*0;
    end 

    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1800 1000]);
     
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', numel(luminancesExamined), ...
           'colsNum', numel(osTimeSteps)+1, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.01);
      
    
    for lumIndex = 1: numel(luminancesExamined) 
        for osTimeStepIndex = 0:numel(osTimeSteps)
            if (osTimeStepIndex > 0)
                theConeMosaic.os.timeStep = osTimeSteps(osTimeStepIndex);
            end
            [isomerizations, photocurrents] = ...
                theConeMosaic.computeForOISequence(theOIsequenceArray{lumIndex}, ...
                'emPaths', emPaths, ...
                'currentFlag', true, ...
                'newNoise', true);

            timeAxis = theConeMosaic.timeAxis + theOIsequenceArray{lumIndex}.timeAxis(1);
     
            % Plot the isomerization signals
            if (osTimeStepIndex == 0)
                subplot('Position', subplotPosVectors(lumIndex, 1).v);
                plot(timeAxis, squeeze(isomerizations(1,1,1,:)), 'r-', 'LineWidth', 1.5);
                hold on;
                plot(timeAxis, squeeze(isomerizations(1,1,2,:)), 'g-', 'LineWidth', 1.5);
                plot(timeAxis, squeeze(isomerizations(1,1,3,:)), 'b-', 'LineWidth', 1.5);
                set(gca, 'XLim', [timeAxis(1) timeAxis(end)]);
                title(sprintf('mean lum: %d cd/m2', luminancesExamined(lumIndex)));
                ylabel(sprintf('absorptions / %d ms', theConeMosaic.integrationTime*1000));
                set(gca, 'YLim', [0 5500]);
                
                if (lumIndex <  numel(luminancesExamined))
                    set(gca, 'XTickLabel', {});
                end
            % Plot the photocurrent signals
            else
                subplot('Position', subplotPosVectors(lumIndex, osTimeStepIndex+1).v);
                plot(timeAxis, squeeze(photocurrents(1,1,1,:)), 'r-', 'LineWidth', 1.5);
                hold on;
                plot(timeAxis, squeeze(photocurrents(1,1,2,:)), 'g-', 'LineWidth', 1.5);
                plot(timeAxis, squeeze(photocurrents(1,1,3,:)), 'b-', 'LineWidth', 1.5);
                set(gca, 'XLim', [timeAxis(1) timeAxis(end)]);
                title(sprintf('os.timeStep = %2.3f ms', theConeMosaic.os.timeStep*1000));
                if (osTimeStepIndex == 1)
                    ylabel('current (pAmps)');
                end
                set(gca, 'YLim', [-90 10]);
                
                if (osTimeStepIndex > 1)
                    set(gca, 'YTickLabel', {});
                end
                
                if (lumIndex == numel(luminancesExamined))
                    xlabel('time (seconds)');
                else
                    set(gca, 'XTickLabel', {});
                end
                set(gca, 'FontSize', 10);
            end        
            grid on; box on; drawnow;
            
            if (osTimeStepIndex == 1)
                % Compute modulation of photocurrents at different mean luminance levels
                for coneIndex = 1:3
                    photocurrent = squeeze(photocurrents(1,1,coneIndex,:));
                    baselineP = photocurrent(end);
                    if (baselineP < 0)
                        deltaPhotocurrent = max(photocurrent - baselineP);
                    else
                        deltaPhotocurrent = 0;
                    end
                    [lumIndex baselineP deltaPhotocurrent] 
                    modulation(lumIndex,coneIndex) = deltaPhotocurrent;
                end % coneIndex
            end
                
        end % osTimeStepIndex
        
    end % lumIndex
    
    % Plot the dynamic range as a function of mean luminance
    figure(2); clf;
    plot(luminancesExamined, modulation(:,1), 'rs-', 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 2.0, 'MarkerSize', 12); hold on
    plot(luminancesExamined, modulation(:,2), 'gs-', 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 2.0, 'MarkerSize', 12);
    plot(luminancesExamined, modulation(:,3), 'bs-', 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 2.0, 'MarkerSize', 12);
    grid on; box on;
    xlabel('luminance (cd/m2)');
    ylabel('modulation');
    set(gca, 'XLim', [70 7000], 'XScale', 'log');
end


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