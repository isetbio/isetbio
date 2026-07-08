function t_fixationalEyeMovementsAndConeMosaicVideo
% Create a hexagonal cone mosaic and fixational eye movement video
%
% Syntax:
%   t_fixationalEyeMovementsAndConeMosaicVideo
%
% Description:
%    Shows how to generate a video of a cone mosaic and a fixational eye
%    movement. Useful for short talks and demonstrations.
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
%    06/21/26  BW   Shorten the video and use the modernized video method.

%
% SkipFile
%

%%
ieInit;

%% Here we use a short integration time to obtain a smoother video.
integrationTime =  0.5 / 1000;  % One half a millisecond (500 usec)

theConeMosaic = cMosaic(...
    'sourceLatticeSizeDegs', 60, ...
    'eccentricityDegs', [5 3], ...
    'sizeDegs', [2 2], ...
    'whichEye', 'right eye', ...
    'eccVaryingConeBlur', true, ...
    'rodIntrusionAdjustedConeAperture', true, ...
    'randomSeed', 1234);

theConeMosaic.integrationTime = integrationTime;

%% Instantiate a fixational eye movement object
% for generating fixational eye movements that include drift and
% microsaccades.
fixEMobj = fixationalEM();

% Generate microsaccades with a mean interval of  150 milliseconds
% Much more often than the default, just for video purposes.
fixEMobj.microSaccadeMeanIntervalSeconds = 0.150;

% Compute one short eye-movement path. At 0.5 msec/sample, a 50 msec trial
% contains 100 samples.
% Here we are fixing the random seed so as to reproduce identical eye
% movements whenever this script is run.
nTrials = 1;
trialLengthSecs = 0.05;
eyeMovementsPerTrial = trialLengthSecs / theConeMosaic.integrationTime;
fixEMobj.computeForCmosaic(theConeMosaic, eyeMovementsPerTrial, ...
    'nTrials', nTrials, 'rSeed', 857);

%% Render the video
% Capture every fourth sample and always capture the final sample. This
% produces 26 frames, or about 0.9 seconds of video at 30 frames/second.
videoExportDir = fullfile(isetbioRootPath,'local');
fixationalEM.generateEMandMosaicComboVideo(...
    fixEMobj, theConeMosaic, ...
    'visualizedFOVdegs', 1.0, ...
    'showMovingMosaicOnSeparateSubFig', true, ...
    'displaycrosshairs', true, ...
    'frameStep', 4, ...
    'frameRate', 30, ...
    'videoFileName', fullfile(videoExportDir, 'eyeMovementVideo.mp4'));

end
