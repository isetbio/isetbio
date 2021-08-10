% Demo using off-axis optics from different data sets
%
% Description:
% Demonstrate usage of @cMosaic, +PolansOptics and +ArtalOptics to generate  
%    eccentricity-varying optics from two different data sets 
%    namely Jaeken and Artal, 2012 (measurements along the horizontal
%    meridian only) and Polans et al, 2015 (measurements along both
%    horizontal and vertical meridians).
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
mosaicSizeDegs = [1 1]*0.15;

%% Mosaic eccentricity (in degrees)
mosaicEccDegsX  = -23:1:22;
mosaicEccDegsX = setdiff(mosaicEccDegsX, [13 14 15 16 17 18]);
mosaicEccDegsY  = 0;

%% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';   

%% choose between {'Polans2015', and 'Artal2012'}
opticsZernikeCoefficientsDataBase = 'Polans2015';            

%% Select ranking of displayed subject
subjectRankOrder = 1;


%% Setup figure
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 2100 1100], 'Color', [1 1 1]);
rowsNum = 4;
colsNum = 10;
sv = NicePlot.getSubPlotPosVectors(...
       'colsNum', colsNum, ...
       'rowsNum', rowsNum, ...
       'heightMargin',  0.08, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.02);
   
%% Generate eccentricity mesh
[X,Y] = meshgrid(mosaicEccDegsX, mosaicEccDegsY);
X = X(:);
Y = Y(:);
R = sqrt(X.^2+Y.^2);

for iEcc = 1:numel(R)
    
    % 
    switch (opticsZernikeCoefficientsDataBase)
        case 'Polans2015'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = PolansOptics.constants.subjectRanking;
            testSubjectID = rankedSujectIDs(subjectRankOrder);

            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

        case 'Artal2012'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
            testSubjectID = rankedSujectIDs(subjectRankOrder);

            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
    end
    
    % Mosaic size: 20 cones across
    conesAcrossMosaic = 20;
    mosaicEccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs([X(iEcc) Y(iEcc)]);
    coneSpacingDegs = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(whichEye, mosaicEccMicrons, 'cones', false);
    sizeDegs = coneSpacingDegs*conesAcrossMosaic*[1 1];
    
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'whichEye', whichEye, ...         
        'sizeDegs', sizeDegs, ...    
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
        'domainVisualizationLimits', domainVisualizationLims, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'labelCones', false, ...
        'noYLabel', true, ...
        'noXlabel', true, ...
        'plotTitle', sprintf('ecc: %d,%d degs\n(spatial support: %2.1f arc min.)', X(iEcc), Y(iEcc), visualizedFraction*sizeDegs(1)*60), ...
        'fontSize', 12);
    hold(ax, 'on');
    cmap = brewermap(1024,'reds');
    alpha = 0.5;
    contourLineColor = [0.0 0.0 0.0];
    cMosaic.semiTransparentContourPlot(ax, psfSupportMicronsX, psfSupportMicronsY, psf, 0.05:0.2:0.95, cmap, alpha, contourLineColor);
    axis(ax, 'square');
end   

NicePlot.exportFigToPDF(sprintf('%s_subject%d.pdf',opticsZernikeCoefficientsDataBase, testSubjectID), hFig, 300);