function t_cMosaicAndOpticsGrid()
% Plot Polans optics across a grid of (x,y) eccentricities
%
% BW: This is a good tutorial to go through with Nicolas and/or David. It
% is a very complicated read of something that should be much simpler. It
% does not use the wvf* methods, but reimplements them in its own way, and
% only their only use is in this context.  I believe that this relies on a
% very important method, oiEnsembleGenerate, that should be at the heart of
% the re-write. We would like the functionality of this tutorial. But we would like it
% to read simply and to rely on the existing tools.
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

%% Control parameters below
%
% Setting fastParamters to true does
% a more limited set of calculations
% and the whole thing runs more quickly.
fastParameters = true;

%% Control saving of figures.
%
% We don't want tutorials saving things into the isetbio source tree
% willy-nilly
saveFigures = false;
figureDir = fullfile(isetbioRootPath,'local',mfilename);
if (saveFigures)
    if (~exist(figureDir,'dir'))
        mkdir(figureDir);
    end
    fprintf('Will save figures into %s\n',figureDir)
else
    fprintf('Not saving figures. Set saveFigures to true in the source to save\n');
end

%% Mosaic size (in degrees)
mosaicSizeDegs = [1 1];

%% Mosaic eccentricity (in degrees)
if (fastParameters)
    fprintf('Running limited set of locations. Change fastParameters to false in the source to get more\n');
    mosaicEccDegsX  = [-12 0 12];
    mosaicEccDegsY  = [-8 0 8];
else
    mosaicEccDegsX  = [-12 -8 -4 -2 0 2 4 8 12];
    mosaicEccDegsY  = [-8 -4 -2 0 2 4 8];
end

%% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';

%% PSF in-focus wavelength
inFocusWavelength = 550;

%% choose between {'Polans2015', and 'Artal2012'}
opticsZernikeCoefficientsDataBase = 'Polans2015';

%% Select ranking of displayed subject and loop
%
% Can provide a list here to see more than one subject.
% For example, use 1:10
for subjectRankOrder = 1

    % Setup a very large and complex figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1600 1200], 'Color', [0 0 0]);

    rowsNum = numel(mosaicEccDegsY);
    colsNum = numel(mosaicEccDegsX);
    sv = NicePlot.getSubPlotPosVectors(...
        'colsNum', colsNum, ...
        'rowsNum', rowsNum, ...
        'heightMargin',   0.06, ...
        'widthMargin',    0.02, ...
        'leftMargin',     0.01, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.05, ...
        'topMargin',      0.02);

    % Generate eccentricity mesh
    [Y,X] = meshgrid(mosaicEccDegsY, mosaicEccDegsX);
    X = X(:);
    Y = Y(:);
    R = sqrt(X.^2+Y.^2);

    % Loop over eccentricities
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
        % sizeDegs = coneSpacingDegs*conesAcrossMosaic*[1 1];

        % Generate mosaic centered at target eccentricity
        cm = cMosaic(...
            'whichEye', whichEye, ...
            'sizeDegs', mosaicSizeDegs, ...
            'eccentricityDegs', [X(iEcc) Y(iEcc)], ...
            'opticalImagePositionDegs', 'mosaic-centered' ...
            );

        % Original comment was brief and not quite right.
        %
        % This call is used to get a PSF.
        % Generate optics appropriate for a particular subject from a
        % particular database at a particular eccentricity and pupil diameter and ...
        [oiEnsemble, psfEnsemble] = ...
            cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
            'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
            'subjectID', testSubjectID, ...
            'pupilDiameterMM', 3.0, ...
            'inFocusWavelength', inFocusWavelength, ...
            'zeroCenterPSF', false, ...
            'flipPSFUpsideDown', false, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'wavefrontSpatialSamples', 501);
        thePSFData = psfEnsemble{1};

        % Visualize PSF
        targetWavelength = 550;
        [~,idx] = min(abs(thePSFData.supportWavelength-targetWavelength));
        psf = squeeze(thePSFData.data(:,:,idx));
        psf = psf/max(psf(:));
        psfSupportMicrons = cm.micronsPerDegree * thePSFData.supportX/60;


        r = floor((iEcc-1)/colsNum);
        r = rowsNum - mod(r,rowsNum);
        c = mod(iEcc-1,colsNum)+1;
        ax = subplot('Position', sv(r,c).v);

        % Ticks and visualization limits
        domainUnits = 'microns';
        
        % Visualize part of the mosaic
        halfWidthMicrons = 20;
        domainVisualizationLims(1:2) = cm.eccentricityMicrons(1) + (halfWidthMicrons+0.5)*[-1 1];
        domainVisualizationLims(3:4) = cm.eccentricityMicrons(2) + (halfWidthMicrons+0.5)*[-1 1];
        domainVisualizationTicks = struct(...
            'x',  sign(cm.eccentricityMicrons(1)) * round(abs(cm.eccentricityMicrons(1))) + halfWidthMicrons*[-1 0 1], ...
            'y',  sign(cm.eccentricityMicrons(2)) * round(abs(cm.eccentricityMicrons(2))) + halfWidthMicrons*[-1 0 1]);

        thePSFdataStruct = struct(...
            'supportXmicrons', psfSupportMicrons, ...
            'supportYmicrons', psfSupportMicrons, ...
            'data', psf ...
            );

        cm.visualize('figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', domainUnits, ...
            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'visualizedConeApertureThetaSamples', 32, ...
            'domainVisualizationLimits', domainVisualizationLims, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'labelConesWithIndices', nan, ...
            'conesAlpha', 0.8, ...
            'conesEdgeAlpha', 0.9, ...
            'withSuperimposedPsf', thePSFdataStruct, ...
            'activationColorMap', brewermap(1024, 'blues'), ...
            'noYLabel', (c > 1), ...
            'noXlabel', (r < rowsNum), ...
            'backgroundColor', 0.85*[1 1 1], ...
            'clearAxesBeforeDrawing', false, ...
            'plotTitle', sprintf('%d, %d', X(iEcc), Y(iEcc)), ...
            'plotTitleColor', [0.5 0.5 0.5], ...
            'fontSize', 16, ...
            'verbose', false, ...
            'plotTitle', sprintf('focus: %2.0fnm, target: %2.0fnm', inFocusWavelength, targetWavelength));
    end

    %% Save figure if desired
    if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(figureDir,sprintf('%s_subject%d_rank%d.pdf',opticsZernikeCoefficientsDataBase, testSubjectID, subjectRankOrder)), hFig, 300);
    end

end

end
