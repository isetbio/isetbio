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
%   t_cMosaicGenerate
%   t_cMosaicSingleEyeMovementPath
%   t_cMosaicMultipleEyeMovementPaths
%   t_cMosaicEccDependentAbsorptionEfficacy
%   t_cMosaicBenchMark
%   t_cMosaicOffAxis
%   t_cMosaicFromConeMosaicHex

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate the ring rays stimulus
scene = sceneCreate('rings rays');
scene = sceneSet(scene, 'fov', 1);

%% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

%% Generate the mosaic
cm = cMosaic(...
    'sizeDegs', [1.0 1.0], ...    % SIZE: 1.0 degs (x) 0.5 degs (y)
    'positionDegs', [1 0], ...  % ECC: (0,0)
    'eccVaryingConeBlur', true ...
    );

% fname = fullfile(isetRootPath,'local','testmosaic');
% overwrite = false;
% ofile = cm.save(fname,overwrite)


%
% To get an OI appropriate for this location, use
% 
%     cm.oiEnsembleGenerate


%% Visualize the mosaic
cm.visualize();

%% Compute multiple noisy response instances of cone excitation response

instancesNum = 2;
[allE, allNoisyE] = cm.compute(oi, 'nTrials', instancesNum);

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

%% First, print out some examples using the 'help' method

regionOfInterest('help');

%% Now choose a small circular ROI

roiCircle = regionOfInterest('shape','ellipse',...
    'center',[1 0],...
    'majorAxisDiameter',0.2,...
    'minorAxisDiameter',0.2);

allE = cm.compute(oi);
% [~,~,allExcitations] = cm.excitations('oi',oi);

roiE = cm.excitations('roi',roiCircle,'all excitations',allE);
mean(roiE)

roiLE = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','L');
mean(roiLE)

cm.get('excitations',varargin)
cm.get('excitations','roi',roiCircle)


% Fewer and lower excitations of S-cones.
[roiSE, roiSIdx] = cm.excitations('roi',roiCircle,...
    'all excitations',allE, ...
    'cone type','S');
mean(roiSE)

[roiE, roiIdx, allE] = cm.excitations('roi',roiCircle,'visualize',true,'all excitations',allE);

%% Illustrate some plots

allE = cm.compute(oi);
params = cm.visualize('params');
params.activation = allE;
cm.visualize(params);

cm.plot('horizontal line',allE, 'y deg',0.3);
cm.plot('horizontal line',allE, 'y deg',0.0);
cm.plot('horizontal line',allE, 'y deg',0.3,'thickness',0.05);

% Make a specific ROI
roiLine = regionOfInterest('shape', 'line', ...
    'from', [.5 .2], 'to', [1.5,.2], ...
    'thickness', 0.1);
cm.plot('horizontal line',allE, 'roi',roiLine);

% Just the M cones
roiLine = regionOfInterest('shape', 'line', ...
    'from', [.5 .3], 'to', [1.5,.3], ...
    'thickness', 0.1);
cm.plot('horizontal line',allE, 'cone type','m','roi',roiLine);

%%




%% Or a line.  Maybe we should always use a rect for a line?  But then
% We would not be able to change the orientation easily.
% Do we want to be able to rotate a rect and an ellipse?


roiLine = regionOfInterest('shape','line',...
    'from',[-2 0],...
    'to',[ 2,0]);
   
cm.visualize();

% Would this be good?
% cm.plot('excitations','line',roiLine,'cone type','L');

idx = roiLine.indicesOfPointsInside(cm.coneRFpositionsDegs);
in = ismember(idx,cm.sConeIndices);
idx = idx(in);
excitations = noiseFreeExcitationResponse(idx);
pos = cm.coneRFpositionsDegs(idx);
ieNewGraphWin;
plot(pos,squeeze(excitations),'-ro');
xlabel('Position (deg)'); ylabel('Excitations');
grid on;





%% Visualize responses
hFig = figure(); clf;
set(hFig, 'Position', [100 300 1500 650]);
activationRange = prctile(noisyExcitationResponseInstances(:), [1 99]);

% Noise-free response
ax = subplot(1,2,1);
cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'activation', noiseFreeExcitationResponse, ...
             'activationRange', activationRange, ...
             'verticalActivationSliceEccentricity', -0.2, ...
             'verticalActivationColorBar', true, ...
             'plotTitle',  'noise-free response');

% Loop over noisy response instances
ax = subplot(1,2,2);
for k = 1:instancesNum
    cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
                 'activation', noisyExcitationResponseInstances(k,:,:), ...
                 'activationRange', activationRange, ...
                 'verticalActivationSliceEccentricity', -0.2, ...
                 'verticalActivationColorBar', true, ...
                 'plotTitle', sprintf('noisy response instance (#%d)', k));
end

%%
