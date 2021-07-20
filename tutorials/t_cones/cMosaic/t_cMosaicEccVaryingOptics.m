% Demo using off-axis optics from different data sets
%
% Description:
% Demo using off-axis optics from two different data sets with off-axis
%    measurements, namely Jaeken and Artal, 2012 (measurements along the horizontal
%    meridian only) and Polans et al, 2015 (measurements along both
%    horizontal and vertical meridians).
%
% See Also:
%   t_cMosaicOffAxisDistortion

% History:
%    07/20/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

%% Initialize
ieInit;
clear;
close all;

%% Mosaic size
mosaicSizeDegs = [1 1]*0.15;
%% Mosaic eccentricity
mosaicEccDegs  = [0 0];
%% On the right eye
whichEye = 'right eye';

%% Generate mosaic centered at target eccentricity
cm = cMosaic(...
        'whichEye', whichEye, ...          % choose {from 'right eye', 'left eye'}
        'sizeDegs', mosaicSizeDegs, ...    % SIZE in degs
        'eccentricityDegs', mosaicEccDegs, ...  % ECC in degs
        'opticalImagePositionDegs', 'mosaic-centered' ...
        );
    
%% Generate optics appropriate for the mosaic's eccentricity  
opticsZernikeCoefficientsDataBase = 'Artal2012';             % choose between {'Polans2015', and 'Artal2012'}

hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1500 800], 'Color', [1 1 1]);
rowsNum = 2;
colsNum = 5;
sv = NicePlot.getSubPlotPosVectors(...
       'colsNum', colsNum, ...
       'rowsNum', rowsNum, ...
       'heightMargin',  0.04, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01); 
       
% Generate and plot foveal PSF for ranked subjects
 switch (opticsZernikeCoefficientsDataBase)
     case 'Polans2015'
         subjectRankOrderList = 1:10;
     case 'Artal2012'
         subjectRankOrderList = 1:5:50;
 end
 
 
for subjectRankOrderIndex = 1:numel(subjectRankOrderList)

    subjectRankOrder = subjectRankOrderList(subjectRankOrderIndex);
    
    switch (opticsZernikeCoefficientsDataBase)
        case 'Polans2015'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = PolansOptics.constants.subjectRanking;
            testSubjectID = rankedSujectIDs(subjectRankOrder);

            % Determine if we need to subtract the subject's central refraction to
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

        case 'Artal2012'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
            testSubjectID = rankedSujectIDs(subjectRankOrder);

            % Determine if we need to subtract the subject's central refraction to
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
    end


    % Generate optics appropriate for the mosaic
    [oiEnsemble, psfEnsemble] = ...
                    cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                    'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', 3.0, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', 501);
    thePSFData = psfEnsemble{1};       

    % Visualize PSF
    targetWavelength = 550;
    [~,idx] = min(abs(thePSFData.supportWavelength-targetWavelength));
    psf = squeeze(thePSFData.data(:,:,idx));
    psf = psf/max(psf(:));
    psfSupportMicrons = cm.micronsPerDegree * thePSFData.supportX/60;
    domainUnits = 'microns';
    domainVisualizationTicks = struct('x',  [-20:10:20], 'y', []);
     
    r = floor((subjectRankOrderIndex-1)/colsNum);
    r = mod(r,rowsNum)+1;
    c = mod(subjectRankOrderIndex-1,colsNum)+1;
    ax = subplot('Position', sv(r,c).v);
        
    cm.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', domainUnits, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'labelCones', false, ...
        'noYLabel', true, ...
        'plotTitle', sprintf('Subj.ID: %d (rank:%d)', testSubjectID, subjectRankOrder));
    hold(ax, 'on');
    cmap = brewermap(1024,'reds');
    alpha = 0.5;
    contourLineColor = [0.3 0.3 0.3];
    cMosaic.semiTransparentContourPlot(ax, psfSupportMicrons, psfSupportMicrons, psf, 0.05:0.1:0.95, cmap, alpha, contourLineColor);
end

NicePlot.exportFigToPDF(sprintf('%s_rankedPSFs.pdf',opticsZernikeCoefficientsDataBase), hFig, 300);