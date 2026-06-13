%% Introduction to the @cMosaic (cone mosaic) object
% We illustrate basic methods for the cone mosaic class, @cMosaic.  This class
% creates a cone mosaic arrays and computes the excitations from different optical
% images (OI).
%
% We illustrate a near fovea @cMosaic object and  compute a number of noisy
% response instances to a static stimulus with no fixational eye movements, so
% the computed responses  consist of a single time point.
%
% We also illustrate how to visualize the mosaic and to plot the computed cone
% excitations.
%
% At the end we illustrate how get all the parameters from the @cMosaic  constructor
% so you can reliably create the same @cMosaic.  Or, you can create a second @cMosaic
% with the same parameters and use it to measure the response to a different stimulus.
% See Also:
%%
% * t_cMosaicPlots
% * t_cMosaicGenerate
% * t_cMosaicSingleEyeMovementPath
% * t_cMosaicMultipleEyeMovementPaths
% * t_cMosaicEccDependentAbsorptionEfficacy
% * t_cMosaicBenchMark
% * t_cMosaicOffAxis
% * t_cMosaicFromConeMosaicHex
%%
% For the simpler, rectangular mosaic, see tutorials/cmosaicrect.
%% Generate the ring rays stimulus

ieInit;

scene = sceneCreate('rings rays', 32, 512);
scene = sceneSet(scene, 'fov', 1.5);
%% Generate the mosaic

positionDegs = [1 0];
sizeDegs = [0.5 0.5];
cm = cMosaic(...
    'size degs', sizeDegs, ...
    'position degs', positionDegs, ...
    'ecc varying cone blur', true ...
    );
%%
% We have a small library of pre-computed mosaics.  You can see which sizes
% and positions we have by typing
%%
%
%   mosaicLoad('help')
%
%%
% You can load a stored mosaic and thus avoid the computation this way

cm = mosaicLoad(sizeDegs,positionDegs);
%%
% This is the standard human lens model.  It is a good approximation over the
% central 10 deg (see Visual Encoding:  Principles and Software, in the Circadian
% Rhythms book, Figure 6).

oi = oiCreate('human');
oi = oiCompute(oi,scene,'pad value','mean');
%%
%
%   % If you want to try a specific OI for this visual field position, you can use
%   zCoeffDatabase   = 'Artal2012';
%   eyeside          = 'right';
%   pupilDiamMM      = 3.0;
%   centerpsf        = true;
%   wave             = 400:10:700;
%   ecc              = 1;  %  In degrees
%   subjectRank      = 15; % [10, 40];
%
%   oi = ...
%       oiPosition(zCoeffDatabase, ...
%       'position',positionDegs, ...
%       'pupil diameter', pupilDiamMM, ...
%       'subject rank', subjectRank, ...
%       'wave',wave, ...
%       'eye side', eyeside,...
% 'center psf',centerpsf);
%
%
%% Visualize of the mosaic and the data

cm.plot('mosaic');
%% To compute a single noise free instance of the excitations

allE = cm.compute(oi);

cm.plot('excitations',allE);
%% To compute multiple noisy response instances of cone excitation response

instancesNum = 3;
[~, allNoisyE] = cm.compute(oi, 'nTrials', instancesNum);

for ii=1:instancesNum
    cm.plot('excitations',allNoisyE(ii,1,:));
end
%% Make a line @regionOfInterest (ROI) and show it superimposed on the activations

% For help on ROIs, use this
%    regionOfInterest('help');

roiLine = regionOfInterest('shape', 'line', ...
    'from', [.8 0.1], 'to', [1.2,-0.1], ...
    'thickness', 0.05);

% Show the ROI on top of the activations
cm.plot('roi',allE, 'roi',roiLine);

% Now choose a small circular ROI
roiCircle = regionOfInterest('shape','ellipse',...
    'center',[1.1 0],...
    'majorAxisDiameter',0.1,...
    'minorAxisDiameter',0.1);

cm.plot('roi',allE, 'roi',roiCircle);
%% Retrieving the excitations within an ROI from all the excitations

roiE = cm.excitations('roi',roiCircle,'all excitations',allE);
mean(roiE)
%%
% The mean excitations in the L cones

roiLE = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','L');
mean(roiLE)
%%
% Notice how few excitations there are, by comparison, in the S-cones

[roiSE, roiSIdx] = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','S');
mean(roiSE)
%% Get the excitations and visualize in one call

[roiE, roiIdx, allE] = cm.excitations('roi',roiCircle,'visualize',true,'all excitations',allE);
%% Visualize responses
% This plotting interface is under development (May, 2025).

hFig = ieFigure;
activationRange = prctile(allNoisyE(:), [1 99]);

% Noise-free response
cm.plot('excitations', allE, ...
    'figureHandle', hFig, ...
    'axesHandle', gca, ...
    'plot title',  'noise-free response');
%%
%
%
% Loop over noisy response instances

hFig = ieFigure;

for k = 1:instancesNum
    cm.plot('excitations', allNoisyE(k,:,:), ...
        'figureHandle', hFig, 'axesHandle', gca, ...
        'plotTitle', sprintf('noisy response instance (#%d)', k));
    pause(1);
end

%% Making a reproducible cMosaic
% For some calculations you would like to re-generate a repeatable, rather than
% random, cone mosaic.  You can use the randomSeed slot for that.

% Here is the first mosaic.  The second argument is the list of
% parameters used to create the mosaic.
[cm1, cm1P] = cMosaic(...
    'sizeDegs', [0.5 0.5], ...
    'positionDegs', [1 0], ...
    'eccVaryingConeBlur', true, ...
    'randomSeed', 12 ...
    );

% Here is the second one, same random seed, so it matches
[cm2, cm2P] = cMosaic(...
    'sizeDegs', [0.5 0.5], ...
    'positionDegs', [1 0], ...
    'eccVaryingConeBlur', true, ...
    'randomSeed', 12 ...
    );

% The parameters are equal because, well, we sent in the same
% parameters.
isequal(cm1P,cm2P)

% You can see the two mosaics are the same this way:
cm1.plot('mosaic');
cm2.plot('mosaic');
%% To create the mosaic from specifying the parameters, you can do this
% Sometimes you just want to control all the parameters.  So using the parameter
% list from above, we can create a mosaic this way.

[cm3, cm3P] = cMosaic(cm1P);

% Check that the parameters remained equal.
isequal(cm1P,cm3P)

% Have a look to see they are equal.
cm3.plot('mosaic');
%% You can change a parameter this way and see that it is different.

cm3P.randomSeed = 11;
cm4 = cMosaic(cm3P);
cm4.plot('mosaic');
%% You can also get the default set of parameters this way

% This creates a random mosaic.
cmP = cMosaicParams;
[cm, cmP1] = cMosaic(cmP);
cm.plot('mosaic');

% The random seed is returned - if you use the parameters again, you
% will get the same mosaic.
cmP1.randomSeed
%%