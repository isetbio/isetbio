function t_fixationalEyeMovementsAndConeMosaicVideo
% Create a hexagonal cone mosaic and fixational eye movement video
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

% Cone hexagonal mosaic params.
% Here we use a short integration time to obtain a smoother video.
integrationTime =  0.5 / 1000;  % One half a millisecond (500 usec)

% Instantiate a hexagonal cone mosaic, or load a cone mosaic
%tutorialsDir = fullfile(isetbioRootPath, 'tutorials');
%resourcesDir = ...
%    fullfile(tutorialsDir, 'eyemovement', 'resources');
%load(fullfile(resourcesDir, 'cmEccBased2.0deg.mat'), ...
%    'theConeMosaic');

theConeMosaic = cMosaic(...
    'sourceLatticeSizeDegs', 60, ...
    'eccentricityDegs', [5 3], ...
    'sizeDegs', [2 2], ...
    'whichEye', 'right eye', ...
    'eccVaryingConeBlur', true, ...
    'rodIntrusionAdjustedConeAperture', true, ...
    'randomSeed', 1234);

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
trialLengthSecs = 0.10;
eyeMovementsPerTrial = trialLengthSecs / theConeMosaic.integrationTime;
fixEMobj.computeForCmosaic(theConeMosaic, eyeMovementsPerTrial, ...
    'nTrials', nTrials, 'rSeed', 857);

% Render the video
videoExportDir = fullfile(isetbioRootPath,'local');
fixationalEM.generateEMandMosaicComboVideo(...
    fixEMobj, theConeMosaic, ...
    'visualizedFOVdegs', 1.0, ...
    'showMovingMosaicOnSeparateSubFig', true, ...
    'displaycrosshairs', true, ...
    'videoFileName', fullfile(videoExportDir, 'eyeMovementVideo.mp4'));

end
