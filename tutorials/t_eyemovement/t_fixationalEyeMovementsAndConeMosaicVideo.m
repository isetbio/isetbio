function t_fixationalEyeMovementsAndConeMosaicVideo
% Create a cone mosaic and fixational eye movement video
%
% Syntax:
%   t_fixationalEyeMovementsAndConeMosaicVideo
%
% Description:
%    Shows how to generate a video of a cone mosaic and a fixational eye
%    movement Useful for talks, demos
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

% History
%    02/15/18  npc  Wrote it.
%    07/18/18  jnm  Formatting

    close all

    % Change to true to render the video
    renderAndSaveVideo = ~true;
    
    % Cone hexagonal mosaic params. 
    % Here we use a short integration time to obtain a smoother video.
    integrationTime = 0.25 / 1000;
    
    % Load a previously generated mosaic to same time
    loadMosaic = true;

    % Instantiate a hexagonal cone mosaic, or load a cone mosaic
    if (loadMosaic)
        tutorialsDir = strrep(isetRootPath, 'isettools', 'tutorials');
        resourcesDir = ...
            fullfile(tutorialsDir, 't_eyemovement', 'resources');
        load(fullfile(resourcesDir, 'cmEccBased2.0deg.mat'), ...
            'theConeMosaic');
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
    nTrials = 2;
    trialLengthSecs = 0.15;
    eyeMovementsPerTrial = trialLengthSecs / theConeMosaic.integrationTime;
    fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
        'nTrials', nTrials, 'rSeed', 857);

    % Render the video
    if (renderAndSaveVideo)
        videoExportDir = uigetdir('~/', 'VIDEO EXPORT DIRECTORY');
        fixationalEM.generateEMandMosaicComboVideo(...
            fixEMobj, theConeMosaic, ...
            'visualizedFOVdegs', 0.5, ...
            'showMovingMosaicOnSeparateSubFig', true, ...
            'displaycrosshairs', true, ...
            'videoFileName', fullfile(videoExportDir, 'video.mp4'));
    end
end
