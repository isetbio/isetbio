function testCmosaic

    % Mosaic size and ecc
    simulationFOVdegs = [1 1];
    simulationEccDegs = [0.2 0.2];
    simulationEye = 'right eye';
    simulationLengthSeconds = 600/1000;
    
    % Select stimulus spatial frequency
    stimFrequencyCPD = 10;
    
    visualizeMosaic = true;
    conditionExamined = 'outer segment length';
    conditionExamined = 'ecc-varying aperture light collection';
    
    conditionExamined = 'ecc-varying os-length and aperture';
    %conditionExamined = 'ecc-varying aperture blur';
        
    %conditionExamined = 'macular pigment';
    %conditionExamined = 'macular pigment with eye movements';
    
    %conditionExamined = 'multiple computes';
    %conditionExamined = 'benchmark1';
    
    [theConeMosaic, theOI] = generateMosaic(...
        simulationEccDegs, ...
        simulationFOVdegs, ...
        simulationEye, ...
        simulationLengthSeconds, ...
        visualizeMosaic);
    
  
    
    % Ask the mosaic to suggest a desired scene pixel size (in
    % degrees), given the state of the 'eccVaryingConeBlur' flag and the 
    % apertures of the mosaic's cones.
    eccVaryingConeBlur = ~true;
    minPixelSizeDegs = theConeMosaic.suggestedScenePixelSizeDegs(eccVaryingConeBlur);
    fprintf('For this mosaic and state of the ''eccVaryingConeBlur'' flag, the suggested scene pixel size is less than or equal to %2.3f arc min\n', ...
        minPixelSizeDegs * 60);
    
    % Load scene
    d = load(fullfile(isetbioRootPath, sprintf('isettools/ganglioncells/demos/scene%dcpd_025.mat', stimFrequencyCPD)), 'theScene');
    theScene = d.theScene;
    
    
    sceneHorizontalFOV = sceneGet(theScene, 'hfov');
    sceneCols = sceneGet(theScene, 'cols');
    fprintf('The pixel size of the loaded scene is %2.3f arc min\n', sceneHorizontalFOV/sceneCols*60);
    
    % Compute the optical image
    theOI = oiCompute(theScene, theOI);
    
    % Compute responses and visualize results
    switch (conditionExamined)
        case 'benchmark1'
            trialsNum = 1000;
            testBenchmark1(theConeMosaic, theOI, trialsNum, 'new cone mosaic',1);
            
        case 'multiple computes'
            testEccVaryingConeBlur(theConeMosaic, theOI,  sprintf('%s_%dCPD', 'ecc-varying aperture blur', stimFrequencyCPD), 1);
            testEccVaryingConeApertureCollection(theConeMosaic, theOI,  'ecc-varying aperture light collection', 1);
            testEccVaryingOuterSegmentLength(theConeMosaic, theOI,  'outer segment length',1);
            testEccVaryingMacularPigmentDensity(theConeMosaic, theOI, 'macular pigment', 1);
            
        case 'ecc-varying aperture light collection'
            % Single aperture blur (median)
            % Comparing collecting area: median vs ecc-changing
            testEccVaryingConeApertureCollection(theConeMosaic, theOI,  conditionExamined,1);
        
        case 'ecc-varying aperture blur'
            % Aperture collecting area increases continuously.
            % Comparing single aperture blur vs 12 - zones of blur
            testEccVaryingConeBlur(theConeMosaic, theOI,  sprintf('%s_%dCPD', conditionExamined,stimFrequencyCPD), 1);
            
        case 'outer segment length'
            % Single aperture blur and collection area (median)
            % Comparing foveal outer segment length vs ecc-dependent outer segment length
            testEccVaryingOuterSegmentLength(theConeMosaic, theOI,  conditionExamined, 2);
            
        case 'ecc-varying os-length and aperture'
            % Comparing foveal outer segment length/median aperture vs ecc-dependent OS length & aperture
            testEccVaryingOSlengthAndAperture(theConeMosaic, theOI,  sprintf('%s_%dCPD', conditionExamined,stimFrequencyCPD), 1);
            
        case 'macular pigment'
            % Comparing foveal macular pigment density vs ecc-dependent density
            testEccVaryingMacularPigmentDensity(theConeMosaic, theOI, sprintf('%s_%dCPD', conditionExamined,stimFrequencyCPD), 1);
            
        case 'macular pigment with eye movements'
            % Comparing foveal macular pigment density vs ecc-dependent
            % density in the presence of eye movements. This is relevant
            % because the MP correction is done on the optical image (to
            % account for the spectral component), although the MP is a
            % property of the mosaic and the mosaic moves with respect to
            % the optical image
            testEccVaryingMacularPigmentDensityWithEyeMovements(theConeMosaic, theOI,  sprintf('%s_%dCPD', conditionExamined,stimFrequencyCPD), 1);
            
             
        otherwise
            error('Unknown condition: ''%s''.', conditionExamined)
    end
end


function testBenchmark1(theConeMosaic, theOI, trialsNum, outputFileName, runNo)
    % Ecc-dependent light collecting area
    theConeMosaic.eccVaryingConeAperture = true;
    % Ecc-dependent OS length
    theConeMosaic.eccVaryingOuterSegmentLength = true;
    % Ecc-dependent macular pigment density
    theConeMosaic.eccVaryingMacularPigmentDensity = true;
    % No ecc-dependent cone blur
    theConeMosaic.eccVaryingConeBlur = true;
    
    % Compute
    tStart = clock();
    [noiseFreeAbsorptions, absorptionInstances] = theConeMosaic.compute(theOI, 'nTrials', trialsNum);
    fprintf('Tile lapsed to compute %d trials: %2.2f seconds\n', trialsNum, etime(clock, tStart));
    
    
    % Compute an estimate of the noise-free absorptions by averaging all response instances
    averagedAbsorptions = mean(absorptionInstances,1);
    
    diffRange = 'bipolar';
    diffColorMap = brewermap(1024, '*RdBu');
    diffScaling = 'linear';
    visualizeResults(theConeMosaic, noiseFreeAbsorptions,  averagedAbsorptions,  ...
            'noise-free', 'averaged fron noisy', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));

    figure(200);
    clf;
    oiResMicrons = oiGet(theOI, 'height spatial resolution')*1e6;
    oiSize = oiGet(theOI, 'size');
    
    cMap = colormap(brewermap(1024, '*greys'));
    
    for coneZoneIndex = 1:(numel(theConeMosaic.blurApertureDiameterMicronsZones)+1)
        ax = subplot(2,5,coneZoneIndex);
       
        if (coneZoneIndex > numel(theConeMosaic.blurApertureDiameterMicronsZones))
            kernel = cMosaic.generateApertureKernel(median(theConeMosaic.blurApertureDiameterMicronsZones), oiResMicrons);
        else
            kernel = cMosaic.generateApertureKernel(theConeMosaic.blurApertureDiameterMicronsZones(coneZoneIndex), oiResMicrons);
        end
        
        kernel = kernel / max(kernel(:));
        
        xx = 1:size(kernel,2);
        yy = 1:size(kernel,1);
        xx = xx - mean(xx);
        yy = yy - mean(yy);
        
        imagesc(ax,xx,yy,kernel);
        axis(ax, 'image');
        set(ax, 'XTick', -10:10, 'YTick', -10:10, 'XLim', 7*[-1 1], 'YLim', 7*[-1 1]);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {}, 'CLim', [-1 1], 'Color', cMap(512,:));
        grid(ax, 'on');
        if (coneZoneIndex > numel(theConeMosaic.blurApertureDiameterMicronsZones))
            title(ax, sprintf('median cone \naperture: %2.1f um', median(theConeMosaic.blurApertureDiameterMicronsZones)));
        else
            title(ax, sprintf('zone %d\naperture: %2.1f um', coneZoneIndex,  theConeMosaic.blurApertureDiameterMicronsZones(coneZoneIndex)));
        end
    end
    colormap(cMap);
        
end

function testEccVaryingMacularPigmentDensityWithEyeMovements(theConeMosaic, theOI, outputFileName, runNo)

    theConeMosaic.eccVaryingConeAperture = true;
    theConeMosaic.eccVaryingOuterSegmentLength = true;
    
    % Ecc-dependent macular pigment density - account for eye movements
    theConeMosaic.eccVaryingMacularPigmentDensity = true;
    theConeMosaic.eccVaryingMacularPigmentDensityDynamic = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI, ...
        'withFixationalEyeMovements', true);
    
    % Foveal macular pigment density - do not account for eye movements
    theConeMosaic.eccVaryingMacularPigmentDensity = true;
    theConeMosaic.eccVaryingMacularPigmentDensityDynamic  = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI, ...
        'withFixationalEyeMovements', true);
    
    diffRange = 'bipolar'; % [-1 1]; %
    diffColorMap = brewermap(1024, '*RdBu');
    diffScaling = 'linear';
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed,  ...
            'ecc-varying MP density, account for EM', 'ecc-varying MP density, ignore EM', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
end




function testEccVaryingMacularPigmentDensity(theConeMosaic, theOI,  outputFileName, runNo)
    theConeMosaic.eccVaryingConeAperture = true;
    theConeMosaic.eccVaryingOuterSegmentLength = true;
    
    % Ecc-dependent macular pigment density
    theConeMosaic.eccVaryingMacularPigmentDensity = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI);
    
    % Foveal macular pigment density
    theConeMosaic.eccVaryingMacularPigmentDensity = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI);
    
    diffRange = 'bipolar';
    diffRange = [-2 2];
    diffColorMap = brewermap(1024, '*RdBu');
    %diffRange = 'unipolar';
    %diffColorMap = brewermap(1024, '*spectral');
    
    diffScaling = 'linear';
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed,  ...
            'ecc-varying MP density', 'mosaic center MP density', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
        
end

function  testEccVaryingOSlengthAndAperture(theConeMosaic, theOI,  outputFileName, runNo)
    
    % Ecc-varying OS length & aperture
    theConeMosaic.eccVaryingConeBlur = ~true;
    theConeMosaic.eccVaryingConeAperture = true;
    theConeMosaic.eccVaryingOuterSegmentLength = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI);
    
    % Foveal OS length
    theConeMosaic.eccVaryingConeBlur = ~true;
    theConeMosaic.eccVaryingConeAperture = ~true;
    theConeMosaic.eccVaryingOuterSegmentLength = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI);
    
    diffRange = 'bipolar';
    diffScaling = 'linear';
    diffColorMap = brewermap(1024, '*RdBu');
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed,  ...
            'ecc-varying OS length & IS area', 'median OS length & median IS area', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
end


function testEccVaryingOuterSegmentLength(theConeMosaic, theOI,  outputFileName, runNo)
    theConeMosaic.eccVaryingConeBlur = ~true;
    theConeMosaic.eccVaryingConeAperture = ~true;
    
    % Ecc-varying OS length
    theConeMosaic.eccVaryingOuterSegmentLength = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI);
    
    % Foveal OS length
    theConeMosaic.eccVaryingOuterSegmentLength = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI);
    
    diffRange = 'bipolar';
    diffScaling = 'linear';
    diffColorMap = brewermap(1024, '*RdBu');
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed,  ...
            'ecc-varying outer segment length', 'median outer segment length', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
        
end


function testEccVaryingConeApertureCollection(theConeMosaic, theOI,  outputFileName, runNo)
   
    theConeMosaic.eccVaryingConeBlur = ~true;
    theConeMosaic.eccVaryingOuterSegmentLength = ~true;
    theConeMosaic.eccVaryingMacularPigmentDensity = ~true;
    
    % Ecc-varying aperture collection
    theConeMosaic.eccVaryingConeAperture = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI);
                
    % Median aperture collection
    theConeMosaic.eccVaryingConeAperture = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI);
    
    diffRange = 'bipolar';
    diffScaling = 'linear';
    diffColorMap = brewermap(1024, '*RdBu');
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed,  ...
            'ecc-varying aperture collection area', 'median aperture collection area', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
        
end

function testEccVaryingConeBlur(theConeMosaic, theOI,  outputFileName, runNo)

    theConeMosaic.eccVaryingOuterSegmentLength = true;
    theConeMosaic.eccVaryingConeAperture = true;
    
    % Ecc-varying blur
    theConeMosaic.eccVaryingConeBlur = true;
    absorptionsCountEcc = theConeMosaic.compute(theOI);
                
    % Median aperture blur
    theConeMosaic.eccVaryingConeBlur = ~true;
    absorptionsCountFixed = theConeMosaic.compute(theOI);


    diffRange = 'bipolar';
    diffScaling = 'linear';
    diffColorMap = brewermap(1024, '*RdBu');
    visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed, ...
            'ecc-varying blur', 'median aperture blur', diffRange,  diffScaling, diffColorMap, sprintf('%s_runNo%d',outputFileName, runNo));
        
end



%     
%     
%     
%     
%     
%     compareApertureAndOSLengthChange = false;
%     compareOSLengthChange = false;
%     compareApertureChange = false;
%     
%     dynamics = 'single instance, single time point';
%     dynamics = 'full';
%     
%     %dynamics = 'zero em, 2 instances, 4 time points';
%     %dynamics = 'one em path, 3 instances';
%     
%     
%     switch (dynamics)
%         case {'one em path, 3 instances', 'full'}
%             
%             if (strcmp(dynamics, 'full'))
%                 
%                 if (compareApertureAndOSLengthChange)
%                     % Compare effects of aperture and os length
%                     theConeMosaic.eccVaryingConeBlur = ~true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%                     absorptionsCountFixed = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj);
% 
%                     theConeMosaic.eccVaryingConeBlur = true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = true;
%                     absorptionsCountEcc = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj);
%                 else 
%                     % Compare effect of MP density
%                     theConeMosaic.eccVaryingConeAperture = true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = true;
%                     
%                     theConeMosaic.eccVaryingMacularPigmentDensity = ~true;
%                     theConeMosaic.eccVaryingMacularPigmentDensityWithEyeMovements = ~true;
%                     absorptionsCountFixed = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj);
% 
%                     theConeMosaic.eccVaryingMacularPigmentDensity = true;
%                     theConeMosaic.eccVaryingMacularPigmentDensityWithEyeMovements = true;
%                     absorptionsCountEcc = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj);
%                 end
%                 
%                 
%             else
%                 instancesNum = 3;
%                 if (compareApertureAndOSLengthChange)
%                     % Compare effects of aperture and os length
%                     theConeMosaic.eccVaryingConeBlur = ~true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%                     absorptionsCount = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj, ...
%                         'nTrials', instancesNum);
% 
%                     theConeMosaic.eccVaryingConeBlur = true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = true;
%                     absorptionsCountEcc = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj, ...
%                         'nTrials', instancesNum);
%                 else 
%                     % Compare effect of MP density
%                     theConeMosaic.eccVaryingConeBlur = true;
%                     theConeMosaic.eccVaryingOuterSegmentLength = true;
%                     
%                     theConeMosaic.eccVaryingMacularPigmentDensity = ~true;
%                     absorptionsCountNoMP = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj, ...
%                         'nTrials', instancesNum);
% 
%                     theConeMosaic.eccVaryingMacularPigmentDensity = true;
%                     absorptionsCount = theConeMosaic.compute(theOI, ...
%                         'fixationalEyeMovementsOBJ', fixEMobj, ...
%                         'nTrials', instancesNum);
%                 end
%             end
%                
%         case 'zero em, 2 instances, 4 time points'
%             timeAxisLength = 4;
%             instancesNum = 2;
%             % Compute absorption density with blur based on median cone aperture 
%             theConeMosaic.eccVaryingConeBlur = ~true;
%             theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%             absorptionsCount = theConeMosaic.compute(theOI, ...
%                 'nTimePoints', timeAxisLength, ...
%                 'nTrials', instancesNum);
%             
%             % Compute absorption density with blur based on ecc-dependent cone aperture 
%             theConeMosaic.eccVaryingConeBlur = true;
%             theConeMosaic.eccVaryingOuterSegmentLength = true;
%             absorptionsCountEcc = theConeMosaic.compute(theOI, ...
%                 'nTimePoints', timeAxisLength, ...
%                 'nTrials', instancesNum);
%             
% 
%         case 'single instance, single time point'
%             
%             if (compareApertureAndOSLengthChange)
%                 % Include effect of aperture change and OS length change
%                 theConeMosaic.eccVaryingOuterSegmentLength = true;
%                 theConeMosaic.eccVaryingConeAperture = true;
%                 
%                 theConeMosaic.eccVaryingConeBlur = true;
%                 absorptionsCountEcc = theConeMosaic.compute(theOI);
%                 
%                 % No aperture change and no OS length change
%                 theConeMosaic.eccVaryingConeBlur = true;
%                 theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%                 absorptionsCountFixed = theConeMosaic.compute(theOI);
%                 
%             elseif (compareOSLengthChange)
%                 % Include effect of aperture change and OS length change
%                 theConeMosaic.eccVaryingConeBlur = true;
%                 theConeMosaic.eccVaryingOuterSegmentLength = true;
%                 absorptionsCountEcc = theConeMosaic.compute(theOI);
%                 
%                 % No aperture change and no OS length change
%                 theConeMosaic.eccVaryingConeBlur = true;
%                 theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%                 absorptionsCountFixed = theConeMosaic.compute(theOI);
%             
%             elseif (compareApertureChange)
%                 % Include effect of aperture change and OS length change
%                 theConeMosaic.eccVaryingOuterSegmentLength = true;
%                 theConeMosaic.eccVaryingConeAperture = true;
%                 theConeMosaic.eccVaryingConeBlur = true;
%                 absorptionsCountEcc = theConeMosaic.compute(theOI);
%                 
%                 % No aperture change and no OS length change
%                 theConeMosaic.eccVaryingConeAperture = ~true;
%                 theConeMosaic.eccVaryingConeBlur = ~true;
%                 theConeMosaic.eccVaryingOuterSegmentLength = ~true;
%                 absorptionsCountFixed = theConeMosaic.compute(theOI);
%                 
%             else
%                 % Compare effect of MP density
%                 theConeMosaic.eccVaryingConeAperture = true;
%                 theConeMosaic.eccVaryingOuterSegmentLength = true;
%                 
%                 theConeMosaic.eccVaryingMacularPigmentDensity = ~true;
%                 absorptionsCountFixed = theConeMosaic.compute(theOI, ...
%                     'doNotaccountForMPdensityChangesDueToEyeMovements', ~true);
% 
%                 
%                 theConeMosaic.eccVaryingMacularPigmentDensity = true;
%                 absorptionsCountEcc = theConeMosaic.compute(theOI, ...
%                     'doNotaccountForMPdensityChangesDueToEyeMovements', ~true);
% 
%             end
% 
%     end % switch
% 
%     if (compareApertureAndOSLengthChange)
%         videoFileName = sprintf('instancesEccApertureAndOSLength_%dCPD', stimFrequencyCPD);
%         visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed, fixEMobj, ...
%              'ecc-varying filtering', 'median aperture filtering', videoFileName);
%     elseif (compareOSLengthChange)
%         videoFileName = sprintf('instancesEccOSLength_%dCPD', stimFrequencyCPD);
%         visualizeResults(theConeMosaic, absorptionsCountEcc, absorptionsCountFixed,  fixEMobj, ...
%             'ecc-varying osLength', 'foveal osLength', videoFileName);
%     elseif (compareApertureChange)
%         videoFileName = sprintf('instancesEccAperture_%dCPD', stimFrequencyCPD);
%         visualizeResults(theConeMosaic, absorptionsCountEcc,  absorptionsCountFixed, fixEMobj, ...
%             'ecc-varying aperture', 'median aperture', videoFileName);
%     else
%         videoFileName = sprintf('instancesMPfovealVsEcc_%dCPD',stimFrequencyCPD);
%         visualizeResults(theConeMosaic, absorptionsCountEcc, absorptionsCountFixed, fixEMobj, ...
%             'ecc-varying MP', 'foveal MP', videoFileName);
%     end
%     
%            
% end


function visualizeResults(theConeMosaic, absorptionsCountCond1, absorptionsCountCond2,  ...
    cond1, cond2, diffRange, diffScaling, diffColorMap, videoFileName)

    if (ischar(diffRange))
        computeDiffRangeFromData = true;
    else
        computeDiffRangeFromData = false;
    end
    
    nTrials = size(absorptionsCountCond1,1);
    nTimePoints = size(absorptionsCountCond1,2);

    m = prctile(absorptionsCountCond1(:), [1 99]);
    
    
    if (nTimePoints > 1)
        videoOBJ = VideoWriter(sprintf('%s.mp4',videoFileName), 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end
    
    for iTrial = 1:nTrials
        for timePoint = 1:nTimePoints
            hFig = figure(100);
            set(hFig, 'Position', [10 10 1800 660], 'Name', sprintf('trial %d, time point: %d', iTrial, timePoint));
            ax = subplot('Position', [0.04 0.07 0.285 0.9]);
            theConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'activation', squeeze(absorptionsCountCond1(iTrial, timePoint,:)), ...
                'activationRange', m, ...
                'horizontalActivationColorBar', true, ...
                'colorBarTickLabelPostFix', 'R*', ...
                'domain', 'microns', ...
                'displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:timePoint), ...
                'crossHairsOnOpticalImageCenter', true, ...
                'fontSize', 16, ...
                'backgroundColor', [0 0 0]);
            if (nTrials > 1)
                 title(ax, sprintf('c1: %s (trial %d/%d)',cond1, iTrial,nTrials));
            else
                 title(ax, sprintf('c1: %s',cond1));
            end
            
            ax = subplot('Position', [0.37 0.07 0.285 0.9]);
            theConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'activation', squeeze(absorptionsCountCond2(iTrial, timePoint,:)), ...
                'activationRange', m, ...
                'horizontalActivationColorBar', true, ...
                'colorBarTickLabelPostFix', 'R*', ...
                'domain', 'microns', ...
                'displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:timePoint), ...
                'crossHairsOnOpticalImageCenter', true, ...
                'fontSize', 16, ...
                'backgroundColor', [0 0 0]);
            title(ax, sprintf('c2: %s',cond2));

            diffMap = 100.0*squeeze(absorptionsCountCond1(iTrial, timePoint,:)-absorptionsCountCond2(iTrial, timePoint,:)) ./ ...
                      squeeze(absorptionsCountCond1(iTrial, timePoint,:));
                  
            if (strcmp(diffScaling, 'log'))
                diffMap = log10(diffMap);
                if (strcmp(diffRange, 'bipolar'))
                    error('Cannot have ''bipolar'' diffRange with ''log'' diffScaling\n');
                end
                colorBarTickLabelPostFix = '';
            else
                colorBarTickLabelPostFix = '%';
            end
            
            if (computeDiffRangeFromData)
                switch (diffRange)
                    case 'unipolar'
                         diffRangeToUse =  prctile(diffMap(:), [0 100]);
                    case 'bipolar'
                         diffRangeToUse =  max(abs(diffMap(:)))*[-1.0 1.0];
                    otherwise
                        error('Unknown value of ''computeDiffRangeFromData'': ''%s''.', diffRange);
                 end
            else
                diffRangeToUse = diffRange;
            end

            
            ax = subplot('Position', [0.70 0.07 0.285 0.9]);
            theConeMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'activation', diffMap, ...
                'activationRange', diffRangeToUse, ...
                'activationColorMap', diffColorMap, ...
                'horizontalActivationColorBar', true, ...
                'colorBarTickLabelPostFix', colorBarTickLabelPostFix, ...
                'domain', 'microns', ...
                'displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:timePoint), ...
                'crossHairsOnOpticalImageCenter', true, ...
                'fontSize', 16, ...
                'backgroundColor', [0 0 0]);
            
            if (strcmp(diffScaling, 'log'))
                title(ax, 'log(100 x (c1 - c2)/ (c1))');
            else
                title(ax, '100 x (c1 - c2)/ (c1)');
            end

            drawnow;
            if (nTimePoints > 1)
                videoOBJ.writeVideo(getframe(hFig));
                NicePlot.exportFigToPDF(sprintf('%s_%d_%d.pdf', videoFileName, iTrial, timePoint), hFig, 300); 
            end
        end
    end
    
    if (nTimePoints > 1)
        videoOBJ.close();
    else
        NicePlot.exportFigToPDF(sprintf('%s.pdf', videoFileName), hFig, 300); 
    end
    
end


function  [theConeMosaic, theOpticsEnsemble] = generateMosaic(...
        simulationEccDegs, simulationFOVdegs, simulationEye, simulationLengthSeconds, visualizeMosaic)
    
    % Cone mosaic
    theConeMosaic = cMosaic(...
            'eccentricityDegs', simulationEccDegs, ...
            'sizeDegs', simulationFOVdegs, ...
            'whichEye', simulationEye, ...
            'coneDensities', [0.62 0.31 0.07 0.0], ...
            'randomSeed', 1);
 
    % Generate eye movement data
    theConeMosaic.emGenSequence(simulationLengthSeconds, ...
        'microsaccadeType', 'heatmap/fixation based', ...
        'microsaccadeMeanIntervalSeconds', 0.2, ...
        'nTrials', 1, ...
        'randomSeed', 10);

    % Optics sampling along the horizontal meridian
    oiSamplingGridDegs(:,1) = simulationEccDegs(1); % [-0.20 0 0.2 0.5 1.0];
    % Optics sampling along the vertical meridian
    oiSamplingGridDegs(:,2) = simulationEccDegs(2)*ones(size(oiSamplingGridDegs,1),1);
    
    % Generate ensemble of optics appropriate for the coneMosaic
    theOpticsEnsemble = theConeMosaic.oiEnsembleGenerate(oiSamplingGridDegs, ...
        'zernikeDataBase', 'Polans2015', ...
        'subjectID', 6, ...
        'pupilDiameterMM', 3.0);
    
    % For now pick the oi which is closest to the center of the mosaic
    [~, oiIndex] = min(sum((bsxfun(@minus, oiSamplingGridDegs, theConeMosaic.eccentricityDegs)).^2,2));
    theOpticsEnsemble = theOpticsEnsemble{oiIndex};

    if (visualizeMosaic)
        hFig = figure(1); clf;
        set(hFig, 'Position', [0 0 1200 1200], 'Color', [1 1 1]);
        ax = subplot('Position', [0.06 0.06 0.93 0.93]);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'microns', ...
            'crossHairsOnOpticalImageCenter', true, ...
            'fontSize', 20);
        NicePlot.exportFigToPDF('newConeMosaic.pdf', hFig, 300);  
    end
end

