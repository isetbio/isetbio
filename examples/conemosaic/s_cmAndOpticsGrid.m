function s_cmAndOpticsGrid()
% Plot Polans optics across a grid of (x,y) eccentricities
%
% See Also:
%   t_cMosaicOffAxisDistortion
%   t_cMosaicRankedSubjectsOptics
%

%% Initialize
ieInit;

%% Control parameters below
%
% Setting fastParamters to true does a more limited set of calculations and
% the whole thing runs more quickly. If you want the larger calculation,
% change it to false.  But put it back so isetbioExampleTest runs quickly.
fastParameters = true;

%% Control saving of figures.
%
% We don't want tutorials saving things into the isetbio source tree
% willy-nilly.
%
% Set default for save to false in any case, as for this example at least
% saving crushes Matlab on the autorun.
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
if (fastParameters)
    mosaicSizeDegs = [0.35 0.35];
    wavefrontSpatialSamples = 201;
    useParforForMosaicGeneration = false;
    visualizedConeApertureThetaSamples = 16;
else
    mosaicSizeDegs = [1 1];
    wavefrontSpatialSamples = 501;
    useParforForMosaicGeneration = true;
    visualizedConeApertureThetaSamples = 32;
end

%% Mosaic eccentricity (in degrees)
if (fastParameters)
    fprintf('Running limited set of locations. Change fastParameters to false in the source to get more\n');
    mosaicEccDegsX  = [-30 0 30];
    mosaicEccDegsY  = [-10 0 10];
else
    mosaicEccDegsX  = [-12 -8 -4 -2 0 2 4 8 12];
    mosaicEccDegsY  = [-8 -4 -2 0 2 4 8];
end

% Size of psf plot times 0.5
psfPlotHalfWidthMicrons = 40;

%% Eye: choose {from 'right eye', 'left eye'}
whichEye = 'right eye';

%% PSF in-focus and target wavelength
inFocusWavelength = 550;
targetWavelength = 550;

%% Choose between {'Polans2015', and 'Artal2012'}
opticsZernikeCoefficientsDataBase = 'Polans2015';

%% Select ranking of displayed subject and loop
%
% Can provide a list here to see more than one subject.
% For example, use 1:10
for subjectRankOrder = 1

    % Setup a large grid figure, capped to fit on typical screens.
    hFig = ieFigure; clf;
    originalFigureUnits = hFig.Units;
    hFig.Units = 'pixels';
    set(hFig, 'Position', localFigurePosition([1600 960]), 'Color', [0 0 0]);
    hFig.Units = originalFigureUnits;

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
        rankedSubjectIDs = PolansOptics.constants.subjectRanking;
        testSubjectID = rankedSubjectIDs(subjectRankOrder);

        % Determine if we need to subtract the subject's central refraction
        subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

        % Generate mosaic centered at target eccentricity
        cm = cMosaic(...
            'whichEye', whichEye, ...
            'sizeDegs', mosaicSizeDegs, ...
            'eccentricityDegs', [X(iEcc) Y(iEcc)], ...
            'opticalImagePositionDegs', 'mosaic-centered', ...
            'useParfor', useParforForMosaicGeneration ...
            );

        % Original comment was brief and not quite right.
        %
        % This call is used to get a PSF.
        % Generate optics appropriate for a particular subject from a
        % particular database at a particular eccentricity and pupil diameter and ...
        [~, psfEnsemble] = ...
            cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
            'zernikeDataBase', opticsZernikeCoefficientsDataBase, ...
            'subjectID', testSubjectID, ...
            'pupilDiameterMM', 3.0, ...
            'inFocusWavelength', inFocusWavelength, ...
            'zeroCenterPSF', false, ...
            'flipPSFUpsideDown', false, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples);
        thePSFData = psfEnsemble{1};

        % Visualize PSF
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
        domainVisualizationLims(1:2) = cm.eccentricityMicrons(1) + (psfPlotHalfWidthMicrons+0.5)*[-1 1];
        domainVisualizationLims(3:4) = cm.eccentricityMicrons(2) + (psfPlotHalfWidthMicrons+0.5)*[-1 1];
        domainVisualizationTicks = struct(...
            'x',  sign(cm.eccentricityMicrons(1)) * round(abs(cm.eccentricityMicrons(1))) + psfPlotHalfWidthMicrons*[-1 0 1], ...
            'y',  sign(cm.eccentricityMicrons(2)) * round(abs(cm.eccentricityMicrons(2))) + psfPlotHalfWidthMicrons*[-1 0 1]);

        thePSFdataStruct = struct(...
            'supportXmicrons', psfSupportMicrons, ...
            'supportYmicrons', psfSupportMicrons, ...
            'data', psf ...
            );

        cm.visualize('figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', domainUnits, ...
            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'visualizedConeApertureThetaSamples', visualizedConeApertureThetaSamples, ...
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
            'plotTitleColor', [0.5 0.5 0.5], ...
            'fontSize', 8, ...
            'verbose', false, ...
            'plotTitle', { sprintf('focus: %2.0fnm, target: %2.0fnm', inFocusWavelength, targetWavelength) ; sprintf('x ecc degs %d, y ecc degs %d', X(iEcc), Y(iEcc)) }, ...
            'plotTitleFontSize', 10 ...
            );
    end

    %% Save figure if desired
    if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(figureDir,sprintf('%s_subject%d_rank%d_%d_%d.pdf',opticsZernikeCoefficientsDataBase, testSubjectID, subjectRankOrder, inFocusWavelength, targetWavelength)), hFig, 300);
    end

end

end

function figurePosition = localFigurePosition(desiredSizePixels)
screenSize = get(groot, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
figureWidth = min(desiredSizePixels(1), max(900, screenWidth - 120));
figureHeight = min(desiredSizePixels(2), max(650, screenHeight - 160));
left = max(10, round((screenWidth - figureWidth)/2));
bottom = max(40, round((screenHeight - figureHeight)/2));
figurePosition = [left bottom figureWidth figureHeight];
end
