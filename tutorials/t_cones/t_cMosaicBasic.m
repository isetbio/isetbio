% Demo basic usage of the new @cMosaic object
%
% Description:
%    Shows basic usage of the new cone mosaic class, @cMosaic.
%    Here, we generate an on-axis (zero eccentricity) cMosaic object and
%    compute a number of noisy response instances to a static stimulus with
%    no fixational eye movements, so the computed responses
%    consist of a single time point.
%
% See Also:
%   tls_cMosaicPlots
%
%   t_cMosaicGenerate
%   t_cMosaicSingleEyeMovementPath
%   t_cMosaicMultipleEyeMovementPaths
%   t_cMosaicEccDependentAbsorptionEfficacy
%   t_cMosaicBenchMark
%   t_cMosaicOffAxis
%   t_cMosaicFromConeMosaicHex

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% What is stored?
%{
[noiseFreeExcitations, noisyFreeExcitations] = cm.compute(oi);

[~,noisyExcitations] = cm.compute(oi);
noisyExcitations = cm.compute(oi,'noise flag','on');


noiseFreeExcitations = cm.compute(oi,'noise flag','off');
noisyExcitations = cm.addNoise('nSamples',8);  % Check noise flag as part of this


cm.compute(oi,'compute type',cellArray)
possibleState = cm.compute('params')
%}

%% Generate the ring rays stimulus

ieInit;
clear;
close all;

scene = sceneCreate('rings rays');
scene = sceneSet(scene, 'fov', 1);

oi = oiCreate;
oi = oiCompute(scene, oi);

%% Generate the mosaic
cm = cMosaic(...
    'sizeDegs', [1.0 1.0], ...    % SIZE: 1.0 degs (x) 0.5 degs (y)
    'positionDegs', [1 0], ...  % ECC: (0,0)
    'eccVaryingConeBlur', true ...
    );

%{
fname = fullfile(isetRootPath,'local','testmosaic');
overwrite = false;
ofile = cm.save(fname,overwrite);
load(ofile,'cmosaic'); cm = cmosaic; 
%}
    
% To get an OI appropriate for this location, use
% 
%     cm.oiEnsembleGenerate
%

%% Visualize of the mosaic and the data

cm.visualize();

%% To compute a single noise free instance of the excitations
allE = cm.compute(oi);

cm.plot('excitations',allE);

%% To compute multiple noisy response instances of cone excitation response

instancesNum = 3;
[~, allNoisyE] = cm.compute(oi, 'nTrials', instancesNum);

for ii=1:instancesNum
    cm.plot('excitations',allNoisyE(ii,1,:));
end

%% Make a line regionOfInterest (ROI) and show it superimposed the activations

% For help on ROIs, use this
%    regionOfInterest('help');

roiLine = regionOfInterest('shape', 'line', ...
    'from', [.5 0.2], 'to', [1.5,-0.2], ...
    'thickness', 0.1);

% Show the ROI on top of the activations
cm.plot('roi',allE, 'roi',roiLine);

% Now choose a small circular ROI
roiCircle = regionOfInterest('shape','ellipse',...
    'center',[1.3 0],...
    'majorAxisDiameter',0.2,...
    'minorAxisDiameter',0.2);

cm.plot('roi',allE, 'roi',roiCircle);

%% Retrieving the excitations within an ROI from all the excitations

roiE = cm.excitations('roi',roiCircle,'all excitations',allE);
mean(roiE)

roiLE = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','L');
mean(roiLE)

% Fewer and lower excitations of S-cones.
[roiSE, roiSIdx] = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','S');
mean(roiSE)

%% Get the excitationa dn visualize in one call

[roiE, roiIdx, allE] = cm.excitations('roi',roiCircle,'visualize',true,'all excitations',allE);

%% Visualize responses - NC code.  Thinking what to do

hFig = figure(); clf;
set(hFig, 'Position', [100 300 1500 650]);
activationRange = prctile(allNoisyE(:), [1 99]);

% Noise-free response
ax = subplot(1,2,1);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
    'activation', allE, ...
    'activationRange', activationRange, ...
    'verticalActivationSliceEccentricity', -0.2, ...
    'verticalActivationColorBar', true, ...
    'plotTitle',  'noise-free response');

% Loop over noisy response instances
ax = subplot(1,2,2);
for k = 1:instancesNum
    cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', allNoisyE(k,:,:), ...
        'activationRange', activationRange, ...
        'verticalActivationSliceEccentricity', -0.2, ...
        'verticalActivationColorBar', true, ...
        'plotTitle', sprintf('noisy response instance (#%d)', k));
    pause(1);
end

