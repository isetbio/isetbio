function t_fixationalEyeMovementsAndConeMosaicVideo
% Shows how to generate a video of a cone mosaic and a fixational eye movement     
% Useful for talks, demos
%

% History
%   02/15/18  npc  Wrote it.

    close all
    
    % Cone hexagonal mosaic params
    integrationTime = 0.25/1000;
    loadMosaic = true;
    
    % Instantiate a hexagonal cone mosaic, or load a cone mosaic
    if (loadMosaic)
        tutorialsDir = strrep(isetRootPath, 'isettools', 'tutorials');
        resourcesDir = fullfile(tutorialsDir, 't_eyemovement', 'resources');
        load(fullfile(resourcesDir, 'cmEccBased2.0deg.mat'), 'theConeMosaic');
    else
        fovDegs = 0.7;
        resamplingFactor = 7;
        theConeMosaic = coneMosaicHex(resamplingFactor, ...
            'fovDegs', fovDegs);
    end
    theConeMosaic.integrationTime = integrationTime;
    
    
    
    % Instantiate a fixational eye movement object for generating
    % fixational eye movements that include drift and microsaccades.
    fixEMobj = fixationalEM();
    
    % Generate microsaccades with a mean interval of  150 milliseconds
    % Much more often than the default, just for video purposes.
    fixEMobj.microSaccadeMeanIntervalSeconds = 0.150;
    
    % Compute nTrials of emPaths for this mosaic
    % Here we are fixing the random seed so as to reproduce identical eye
    % movements whenever this script is run.
    nTrials = 2; trialLengthSecs = 0.5;
    eyeMovementsPerTrial = trialLengthSecs/theConeMosaic.integrationTime;
    fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
        'nTrials', nTrials, ...
        'rSeed', 857 ...
        );
    
    % Render the video
    fixationalEM.generateEMandMosaicComboVideo(fixEMobj, theConeMosaic, ...
        'visualizedFOVdegs', 0.5, ...
        'showMovingMosaicOnSeparateSubFig', true, ...
        'displaycrosshairs', true, ...
        'videoFileName', 'video.mp4');
    
end
