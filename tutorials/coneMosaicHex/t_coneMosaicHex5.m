%% t_coneMosaicHex5
%
% Shows how to generate hexagonal mosaic with varying density and
% spatial inhomogeneities.
%
% NPC ISETBIO Team, Copyright 2016

%% Initialize
ieInit; clear; close all;
    
rng('default'); rng(219347);

%% Unit test 1: generate a sptially-varying density hex mosaic positioned at (0.5mm, 0.0mm)
% Mosaic Parameters
mosaicParams = struct(...
    'resamplingFactor', 7, ...          
      'varyingDensity', true, ...
          'centerInMM', [0.5 0.0], ...       % 0.5 mm horizontal, 0.0 vertical
                'size', [32 50] ...          % generate from a rectangular mosaic of 32x50 cones
  );
commandwindow
fprintf('\n<strong>Hit enter to create a spatially-varying density hex mosaic positioned at %2.1f mm, %2.1fmm. </strong>', mosaicParams.centerInMM(1), mosaicParams.centerInMM(2));
pause

% Generate the hex grid
theHexMosaic1 = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.varyingDensity, ...
            'center', mosaicParams.centerInMM*1e-3, ....
              'size', mosaicParams.size ...
    ); 
% Print some grid info and visualize it
theHexMosaic1.displayInfo();
theHexMosaic1.visualizeGrid('overlayConeDensityContour', true, 'generateNewFigure', true);

%% Unit test 2: generate a sptially-varying density hex mosaic positioned at (0.0, 0.5)
mosaicParams.centerInMM = [0.0 0.5];
commandwindow
fprintf('\n<strong>Hit enter to create a spatially-varying density hex mosaic positioned at %2.1f mm, %2.1fmm. </strong>', mosaicParams.centerInMM(1), mosaicParams.centerInMM(2));
pause

theHexMosaic2 = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.varyingDensity, ...
            'center', mosaicParams.centerInMM*1e-3, ....
              'size', mosaicParams.size ...
    ); 
% Print some grid info and visualize it
theHexMosaic2.displayInfo();
theHexMosaic2.visualizeGrid('overlayConeDensityContour', true, 'generateNewFigure', true);

%% Unit test 3: generate a sptially-varying density hex mosaic positioned at (0.1, 0.1)
mosaicParams.centerInMM = [0.1 0.1];
commandwindow
fprintf('\n<strong>Hit enter to create a spatially-varying density hex mosaic positioned at %2.1f mm, %2.1fmm. </strong>', mosaicParams.centerInMM(1), mosaicParams.centerInMM(2));
pause

theHexMosaic3 = coneMosaicHex(mosaicParams.resamplingFactor, mosaicParams.varyingDensity, ...
            'center', mosaicParams.centerInMM*1e-3, ....
              'size', mosaicParams.size ...
    ); 
% Print some grid info and visualize it
theHexMosaic3.displayInfo();
theHexMosaic3.visualizeGrid('overlayConeDensityContour', true, 'generateNewFigure', true);


%% Unit test 4: compute and display activation maps to a Gabor stimulus
commandwindow
fprintf('\n<strong>Hit enter to compute and visualize isomerizations maps for the 3 mosaics for an achromatic Gabor scene. </strong>', mosaicParams.centerInMM(1), mosaicParams.centerInMM(2));
pause
% Load acrhomatic Gabor scene
[dirName,~] = fileparts(which(mfilename()));
load(fullfile(dirName,'GaborAchromScene.mat'))
gaborScene = sceneSet(gaborScene,'fov', 2.0);

% Compute the optical image
oi = oiCreate;
oi = oiCompute(gaborScene,oi); 

% Compute isomerizations for the different mosaics and display the results
isomerizationsGabor1 = theHexMosaic1.compute(oi,'currentFlag',false);
isomerizationsGabor2 = theHexMosaic2.compute(oi,'currentFlag',false);
isomerizationsGabor3 = theHexMosaic3.compute(oi,'currentFlag',false);

theHexMosaic1.visualizeActivationMaps(...
     isomerizationsGabor1, ...                                      % the signal matrix
       'mapType', 'modulated hexagons', ...                          % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
    'signalName', 'isomerizations (R*/cone/integration time)', ...   % colormap title (signal name and units)
      'colorMap', jet(1024), ...                                    % colormap to use for displaying activation level
    'figureSize', [1550 950] ...                                     % figure size in pixels
    );

theHexMosaic2.visualizeActivationMaps(...
     isomerizationsGabor2, ...                                      % the signal matrix
       'mapType', 'modulated hexagons', ...                          % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
    'signalName', 'isomerizations (R*/cone/integration time)', ...   % colormap title (signal name and units)
      'colorMap', jet(1024), ...                                    % colormap to use for displaying activation level
    'figureSize', [1550 950] ...                                     % figure size in pixels
    );

theHexMosaic3.visualizeActivationMaps(...
     isomerizationsGabor3, ...                                      % the signal matrix
       'mapType', 'modulated hexagons', ...                          % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
    'signalName', 'isomerizations (R*/cone/integration time)', ...   % colormap title (signal name and units)
      'colorMap', jet(1024), ...                                    % colormap to use for displaying activation level
    'figureSize', [1550 950] ...                                     % figure size in pixels
    );

