%%t_coneMosaicAbsorptionsEccBased  Basic introduction to the cone mosaic object.
%
% Description:
%    Show how to create a cone mosaic object and compute cone isomerizatoins
%    across a set of small eye movements.  This lets you look at the
%    result in the coneMosaic window.
%

%% Build a simple scene and oi (retinal image) for computing
%
% First the scene
sz          = 256;
barSlope    = 2.6;
fieldOfView = 1;
s = sceneCreate('slantedBar', sz, barSlope, fieldOfView);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi,s);

%% Build a default cone mosaic and compute isomerizatoins
saveMosaic = ~true;
mosaicFileName = sprintf('mosaic.mat');
if (saveMosaic)
    cMosaic  = coneMosaicHex(7, ...
            'name', 'test', ...
            'fovDegs', fieldOfView, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', ~true, ...
            'sConeMinDistanceFactor', 2.0, ... 
            'sConeFreeRadiusMicrons', 0.15*300, ...                   
            'spatialDensity', [0 6/10 3/10 1/10], ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...         
            'latticeAdjustmentDelaunayToleranceF', 0.05, ... 
            'maxGridAdjustmentIterations', 20, ...
            'marginF', [] ... 
        );
    save(mosaicFileName, 'cMosaic', '-v7.3');
else
    load(mosaicFileName, 'cMosaic');
end

cMosaic.visualizeGrid(...
        'generateNewFigure', true, ...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true, ...
        'generateNewFigure', true);
    

%% Generate a sequence of 3 eye posistions.
cMosaic.emGenSequence(3);
cMosaic.emPositions = 0 * cMosaic.emPositions;

%% Compute isomerizations for each eye position.
fprintf('Computing\n')
cMosaic.eccBasedConeQuantalEfficiency = true;
absorptions = cMosaic.compute(oi);
activation = squeeze(absorptions(1,:,:,1));
signalRange = [min(activation(:)) max(activation(:))]

figure(100); clf;
subplot(2,2,1);
cMosaic.renderActivationMap(gca, activation, ...
         'signalRange', signalRange, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0])
title('With ecc-based corrections');



% Now compute with eccBasedConeQuantalEfficiency set to true
cMosaic.eccBasedConeQuantalEfficiency = true;

% Visualize grid
cMosaic.visualizeGrid(...
        'generateNewFigure', true, ...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true, ...
        'generateNewFigure', true);
    

%% Compute isomerizations for each eye position.
fprintf('Computing again\n')
cMosaic.emGenSequence(3);
cMosaic.emPositions = 0 * cMosaic.emPositions;


absorptions = cMosaic.compute(oi);
activation = squeeze(absorptions(1,:,:,1));
signalRange2 = [min(activation(:)) max(activation(:))]

figure(100);
subplot(2,2,2);
cMosaic.renderActivationMap(gca, activation, ...
         'signalRange', signalRange, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0]);
title('With ecc-based corrections (again)');


%% Compute isomerizations for each eye position.
fprintf('Computing again\n')
cMosaic.emGenSequence(3);
cMosaic.emPositions = 0 * cMosaic.emPositions;

cMosaic.eccBasedConeQuantalEfficiency = ~true;
absorptions = cMosaic.compute(oi);
activation = squeeze(absorptions(1,:,:,1));
signalRange3 = [min(activation(:)) max(activation(:))]

figure(100);
subplot(2,2,3);
cMosaic.renderActivationMap(gca, activation, ...
         'signalRange', signalRange3, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0]);
title('Without ecc-based corrections');
