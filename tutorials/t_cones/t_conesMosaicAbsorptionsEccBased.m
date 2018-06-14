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
%
% Create the coneMosaic object
%cMosaic = coneMosaic;  
% Set size to show about half the scene.  Speeds things
% up.
%cMosaic.setSizeToFOV(0.5*sceneGet(s,'fov')); 

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

cMosaic.eccBasedConeQuantalEfficiency = ~true;
cMosaic.visualizeGrid(...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true, ...
        'generateNewFigure', true);
    

%% Generate a sequence of 100 eye posistions.
cMosaic.emGenSequence(2);
cMosaic.emPositions = 0 * cMosaic.emPositions;

%% Compute isomerizations for each eye position.
absorptions = cMosaic.compute(oi);
activation = squeeze(absorptions(1,:,:,1));
signalRange = [min(activation(:)) max(activation(:))]

figure(100);
subplot(1,2,1);
cMosaic.renderActivationMap(gca, activation, ...
         'signalRange', signalRange, ...
         'mapType', 'modulated disks', ...
         'colorMap', gray(1024));
set(gca, 'Color', [0 0 0])
%%

