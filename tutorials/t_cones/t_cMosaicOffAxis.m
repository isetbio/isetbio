% Demo computation with an off-axis  @cMosaic object
%
% Description:
%    Shows how to generate and use the new cone mosaic class, @cMosaic.
%    Here, we generate off-axis (located at different eccentricities)
%    cMosaic objects and compute their mean responses to a static stimulus
%    (rigs/rays) with no eye movements. 
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxisDistortion
%   t_cMosaicBenchMark

% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Generate the ring rays stimulus
scene = sceneCreate('ringsrays', 10, 512);
scene = sceneSet(scene, 'fov', 8.0);

%% Set up figures and subfigs
hFig = figure(1);
set(hFig, 'Position', [10 10 1400 1200]);

sv = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', 3, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.03);

%% Generate mosaics and compute responses on a 3x3 grid of eccentricities
deltaDeg = 0.7;
for xOffset = 1:3
for yOffset = 1:3

    % Mosaic eccentricity
    mosaicEcc = [(xOffset-2)*deltaDeg (yOffset-2)*deltaDeg];
    
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', [1 1]*0.7, ...          % SIZE: 0.7 degs (x) 0.7 degs (y)
        'eccentricityDegs', mosaicEcc ...  % ECC: varying
        );

    % Generate standard human optics
    oi = oiCreate('wvf human');
    
    % Compute the optical image of the scene
    oi = oiCompute(scene, oi);
 
    % Compute the noise-free excitation response
    noiseFreeExcitationResponse = cm.compute(oi, 'opticalImagePositionDegs', [0 0]);

    % Visualize mosaic response
    cm.visualize('figureHandle', hFig, ...
        'axesHandle', subplot('Position', sv(end-yOffset+1,xOffset).v), ...
        'domain', 'degrees', ...
        'crossHairsOnMosaicCenter', true, ...
        'domainVisualizationTicks', struct('x', -5:0.1:5, 'y', -5:0.1:5), ...
        'activation', noiseFreeExcitationResponse, ...
        'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));

end
end

