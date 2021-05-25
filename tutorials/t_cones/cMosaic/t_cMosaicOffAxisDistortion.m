% Demo computation with an off-axis @cMosaic and optics
%
% Description:
%    Shows how to generate and use the new cone mosaic class, @cMosaic.
%    Here, we generate off-axis (located at different eccentricities)
%    cMosaic objects and corresponding optics, and compute mosaic responses
%    to a grid stimulus so at to reveal that distortions introduced by the 
%    optics at different eccentricities
%
% See Also:
%   t_cMosaicBasic
%   t_cMosaicOffAxis
%   t_cMosaicBenchMark


% History:
%    03/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.


%% Initialize
ieInit;
clear;
close all;

%% Mosaic size
mosaicSizeDegs = [1 1]*0.7;

%% Generate the linegrid stimulus
scene = sceneCreate('distortiongrid', 512, 100, 'ep');
scene = sceneSet(scene, 'fov', max(mosaicSizeDegs)*1.1);

%% Set up figures and subfigs
hFig1 = figure(1);
set(hFig1, 'Position', [10 10 1400 1200], 'Name', 'optics');

hFig2 = figure(2);
set(hFig2, 'Position', [1000 10 1400 1200], 'Name', 'mosaic responses');

sv = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', 3, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.03);


%% Select a subject to see the effects of different optics
PolansSubject1 = 9;  % subject with vertically-elongated PSFs
PolansSubject2 = 8;  % subject with  horizontally-elongated PSFs
PolansSubject = PolansSubject1;

%% Generate mosaics and compute responses on a 3x3 grid of eccentricities
deltaDeg = 1.0;
for xOffset = 1:3
for yOffset = 1:3

    % Mosaic eccentricity
    mosaicEcc = [(xOffset-1)*deltaDeg (yOffset-1)*deltaDeg];
    
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', mosaicSizeDegs, ...    % SIZE in degs
        'eccentricityDegs', mosaicEcc, ...  % ECC in degs
        'opticalImagePositionDegs', 'mosaic-centered' ...
        );

    % Generate optics appropriate for the mosaic's eccentricity
    [oiEnsemble, psfEnsemble] = cm.oiEnsembleGenerate(mosaicEcc, ...
        'zernikeDataBase', 'Polans2015', ...
        'subjectID', PolansSubject, ...
        'pupilDiameterMM', 3.0, ...
        'subtractCentralRefraction', true);
    
    % Compute the optical image of the scene
    oi = oiCompute(scene, oiEnsemble{1});
 
    % Compute the noise-free excitation response
    noiseFreeExcitationResponse = cm.compute(oi);

    % Visualize optics
    figure(1);
    thePSF = psfEnsemble{1};
    [~, wIdx] = min(abs(thePSF.supportWavelength-550));
    wavePSF = squeeze(thePSF.data(:,:,wIdx));
    zLevels = 0.1:0.1:0.9;
    xyRangeArcMin = 5*[-1 1];
    PolansOptics.renderPSF(subplot('Position', sv(end-yOffset+1,xOffset).v), ...
        thePSF.supportX, thePSF.supportY, wavePSF/max(wavePSF(:)), ...
        xyRangeArcMin, zLevels,  gray(1024), [0 0 0], ...
        'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));
    
    % Visualize mosaic responses
    figure(2);
    % Visualize mosaic response
    cm.visualize('figureHandle', hFig2, ...
        'axesHandle', subplot('Position', sv(end-yOffset+1,xOffset).v), ...
        'domain', 'degrees', ...
        'domainVisualizationTicks', struct('x', -5:0.1:5, 'y', -5:0.1:5), ...
        'activation', noiseFreeExcitationResponse, ...
        'plotTitle',  sprintf('ecc: %2.1f, %2.1f degs', mosaicEcc(1), mosaicEcc(2)));

end
end

