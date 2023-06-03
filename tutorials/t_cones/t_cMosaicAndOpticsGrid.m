function t_cMosaicAndOpticsGrid()
% Plot Polans optics across a grid of (x,y) eccentricities
%
% Description:
%
% See Also:
%   t_cMosaicOffAxisDistortion
%   t_cMosaicRankedSubjectsOptics

% History:
%    07/20/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

%% Initialize
ieInit;
clear;
close all;

%% Mosaic size (in degrees)
mosaicSizeDegs = [1 1]*0.5;

%% Mosaic eccentricity (in degrees)
mosaicEccDegsX  = [-12 -8 -4 -2 0 2 4 8 12];
mosaicEccDegsY  = [-8 -4 -2 0 2 4 8];

%% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';   

%% choose between {'Polans2015', and 'Artal2012'}
opticsZernikeCoefficientsDataBase = 'Polans2015';            

%% Select ranking of displayed subject
for subjectRankOrder = 1:10


%% Setup figure
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1600 1200], 'Color', [0 0 0]);

rowsNum = numel(mosaicEccDegsY);
colsNum = numel(mosaicEccDegsX);
sv = NicePlot.getSubPlotPosVectors(...
       'colsNum', colsNum, ...
       'rowsNum', rowsNum, ...
       'heightMargin',  0.015, ...
       'widthMargin',    0.005, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.01, ...
       'topMargin',      0.01);
   
%% Generate eccentricity mesh
[Y,X] = meshgrid(mosaicEccDegsY, mosaicEccDegsX);
X = X(:);
Y = Y(:);
R = sqrt(X.^2+Y.^2);

for iEcc = 1:numel(R)
    
    % Obtain subject IDs ranking in decreasing foveal resolution
    rankedSujectIDs = PolansOptics.constants.subjectRanking;
    testSubjectID = rankedSujectIDs(subjectRankOrder);

    % Determine if we need to subtract the subject's central refraction
    subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    
    % Mosaic size: 20 cones across
    conesAcrossMosaic = 20;
    mosaicEccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs([X(iEcc) Y(iEcc)]);
    coneSpacingDegs = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(whichEye, mosaicEccMicrons, 'cones', false);
    sizeDegs = coneSpacingDegs*conesAcrossMosaic*[1 1];
    
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'whichEye', whichEye, ...         
        'sizeDegs', mosaicSizeDegs, ...    
        'eccentricityDegs', [X(iEcc) Y(iEcc)], ...  
        'opticalImagePositionDegs', 'mosaic-centered' ...
        );
    
    % Generate optics appropriate for the mosaic's eccentricity  
    [oiEnsemble, psfEnsemble] = ...
                    cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                    'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', 3.0, ...
                    'zeroCenterPSF', false, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', 501);
    thePSFData = psfEnsemble{1};       

    % Visualize PSF
    targetWavelength = 550;
    [~,idx] = min(abs(thePSFData.supportWavelength-targetWavelength));
    psf = squeeze(thePSFData.data(:,:,idx));
    psf = psf/max(psf(:));
    psfSupportMicrons = cm.micronsPerDegree * thePSFData.supportX/60;
    psfSupportMicronsX = psfSupportMicrons + cm.eccentricityMicrons(1);
    psfSupportMicronsY = psfSupportMicrons + cm.eccentricityMicrons(2);
    
    
    r = floor((iEcc-1)/colsNum);
    r = mod(r,rowsNum)+1;
    c = mod(iEcc-1,colsNum)+1;
    ax = subplot('Position', sv(r,c).v);
        
    % Ticks and visualization limits
    domainUnits = 'microns';
    domainVisualizationTicks = struct('x',  [nan], 'y', [nan]);
    % Visualize half of the mosaic
    visualizedFraction = 0.3;
    w = cm.sizeDegs(1)*cm.micronsPerDegree;
    domainVisualizationLims(1:2) = cm.eccentricityMicrons(1) + visualizedFraction*w/2*[-1 1];
    domainVisualizationLims(3:4) = cm.eccentricityMicrons(2) + visualizedFraction*w/2*[-1 1];
    
    cm.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', domainUnits, ...
        'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'domainVisualizationLimits', domainVisualizationLims, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'labelCones', true, ...
        'conesAlpha', 0.0, ...
        'conesEdgeAlpha', 0.4, ...
        'noYLabel', true, ...
        'noXlabel', true, ...
        'backgroundColor', [1 1 1], ...
        'clearAxesBeforeDrawing', false, ...
        'plotTitle', sprintf('%d, %d', X(iEcc), Y(iEcc)), ...
        'plotTitleColor', [0.5 0.5 0.5], ...
        'fontSize', 12);
    hold(ax, 'on');
    cmap = brewermap(1024,'reds');
    alpha = 0.5;
    contourLineColor = [0.0 0.0 0.0];
    cMosaic.semiTransparentContourPlot(ax, psfSupportMicronsX, psfSupportMicronsY, psf, 0.05:0.2:0.95, cmap, alpha, contourLineColor);
    axis(ax, 'square');
end   

NicePlot.exportFigToPDF(sprintf('%s_subject%d_rank%d.pdf',opticsZernikeCoefficientsDataBase, testSubjectID, subjectRankOrder), hFig, 300);
end

end
