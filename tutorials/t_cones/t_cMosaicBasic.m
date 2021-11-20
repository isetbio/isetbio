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
scene = sceneSet(scene, 'fov', 0.5);

%% Compute the optical image
oi = oiCreate;
oi = oiCompute(scene, oi);

%% Generate the mosaic
cm = cMosaic(...
    'sizeDegs', [1.0 1.5]*2, ...    % SIZE: 1.0 degs (x) 0.5 degs (y)
    'eccentricityDegs', [0 0], ...  % ECC: (0,0)
    'eccVaryingConeBlur', true ...
    );

%% Visualize the mosaic
cm.visualize();

%% Compute 8 noisy response instances of cone excitation response
instancesNum = 2;
[noiseFreeExcitationResponse, noisyExcitationResponseInstances] = cm.compute(oi, ...
    'nTrials', instancesNum);

%% First, print out some examples using the 'help' method

regionOfInterest('help');

%% Now choose a small circular ROI

roiCircle = regionOfInterest('shape','ellipse',...
    'center',[0 0],...
    'majorAxisDiameter',0.2,...
    'minorAxisDiameter',0.2);

% Find indices of the cones whose positions are within the ROI
idx = roiCircle.indicesOfPointsInside(cm.coneRFpositionsDegs);

% Find the excitations at those positions
excitations = noiseFreeExcitationResponse(idx);
mean(excitations)

in = ismember(idx,cm.lConeIndices);
idxL = idx(in);
mean(noiseFreeExcitationResponse(idxL))

% M
in = ismember(idx,cm.mConeIndices);
idxM = idx(in);
mean(noiseFreeExcitationResponse(idxM))

% Good, really.  No S-cones.
in = ismember(idx,cm.sConeIndices);
idxS = idx(in);
mean(noiseFreeExcitationResponse(idxS))

% Show the region of the extraction.
testActivation = noiseFreeExcitationResponse;
testActivation(idx) = 1000;

% This is the region the data come from
params = cm.visualize('params');
params.activation = testActivation;
params.verticalActivationColorBar = true;
cm.visualize(params);

params.activation = noiseFreeExcitationResponse;
cm.visualize(params);



%% Or a line

roiLine = regionOfInterest('shape','line',...
    'from',[-0.5 0],...
    'to',[ 0.5,0]);
   
idx = roiLine.indicesOfPointsInside(cm.coneRFpositionsDegs);
excitations = noiseFreeExcitationResponse(idx);



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
