% Script to examine the max Strehl optics of a subject at different eccentricities
% This script is also used to generate materials for the validation
% figures of our PLOS2024 paper
%
% Usage:
%{
    t_testOptimizedStrehlOptics  
%}

% Which optics to use
whichZernikeDataBase = 'Polans2015';
whichSubjectID = 2;
whichEye = 'right eye';
noLCA = false;

% Root directory
intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();
intermediateDataDir = fullfile(intermediateDataDir, 'SLIM', 'demos');
matFile = fullfile(intermediateDataDir, sprintf('StrehlRatioAnalysis_%s_Subj_%d.mat', whichZernikeDataBase,whichSubjectID));
%matFile = fullfile(intermediateDataDir, sprintf('StrehlRatioAnalysis_%s_Subj_%d_HiResGrid.mat', whichZernikeDataBase,whichSubjectID));

ASPPManager = AppleSiliconParPoolManager('half max');

opticsParams = struct(...
    'ZernikeDataBase', whichZernikeDataBase, ...
    'subjectID', whichSubjectID, ...
    'pupilDiameterMM', 3.0, ...
    'noLCA', noLCA, ...
    'zeroCenterPSF', true, ...
    'modification', 'maxStrehlRatio');

% Set cone aperture modifiers
sigmaGaussian = 0.204;
coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');


inFocusWavelength = 550;
inFocusWavelengthIndex = [];

recomputeMaxStrehlRatioPSFs = ~true;

if (recomputeOptimalStrehlRatioPSFs)
    sampledXpos = [-10 -5 -2 -1 0 1 2 5 10]
    sampledYpos = [-5 -1 0 1 5]

    thePSFs = cell(numel(sampledXpos), numel(sampledYpos));
    theConeMosaics = cell(numel(sampledXpos), numel(sampledYpos));
    theOptimalStrehlRatioDefocusDiopters = zeros(numel(sampledXpos), numel(sampledYpos));
    theOptimalStrehlRatio = zeros(numel(sampledXpos), numel(sampledYpos));

    for ix = 1:numel(sampledXpos)
        parfor iy = 1:numel(sampledYpos)
            fprintf('Generating maxStrehl ratio optics at %2.1f,%2.1f\n', sampledXpos(ix), sampledYpos(iy))
            eccentricityPosDegs = [sampledXpos(ix) sampledYpos(iy)];
            theConeMosaics{ix,iy} = cMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'whichEye', whichEye, ...
                'eccentricityDegs', eccentricityPosDegs, ...
                'sizeDegs', [1 1], ...
                'overlappingConeFractionForElimination', 0.01, ...
                'rodIntrusionAdjustedConeAperture', true, ...
                'coneApertureModifiers', coneApertureModifiers ...
                );

            theConeMosaic = theConeMosaics{ix,iy};
            [~, thePSFatCurrentEccentricity, theOptimalStrehlRatioDefocusDiopters(ix,iy), theOptimalStrehlRatio(ix,iy)] = ...
                RGCMosaicConstructor.helper.optics.generate(theConeMosaic, eccentricityPosDegs, opticsParams, ...
                    'visualizeStrehlRatioOptimization', false);

            thePSFs{ix,iy} = thePSFatCurrentEccentricity;
         end
    end
    % Save the Strehl ratio maximized STFs
    fprintf('Saving computed max Strehl ratio PSFs to\n%s\n', matFile);
    save(matFile, 'sampledXpos', 'sampledYpos', 'thePSFs', ...
        'theConeMosaics', 'theOptimalStrehlRatioDefocusDiopters', 'theOptimalStrehlRatio', '-v7.3');
else
    % Load the Strehl ratio maximized STFs
    fprintf('Loading previously computed max Strehl ratio PSFs from\n%s\n', matFile);
    load(matFile, 'sampledXpos', 'sampledYpos', 'thePSFs', ...
        'theConeMosaics', 'theOptimalStrehlRatioDefocusDiopters', 'theOptimalStrehlRatio');
end


wavelengthSupport = theConeMosaics{1,1}.wave;
[achromaticIncrementRadiancePhotons, achromaticDecrementRadiancePhotons] = ...
    RGCMosaicConstructor.helper.simulateExperiment.stimulusDifferentialSpectrum(wavelengthSupport, 'Achromatic');

[LconeIncrementRadiancePhotons, LconeDecrementRadiancePhotons] = ...
    RGCMosaicConstructor.helper.simulateExperiment.stimulusDifferentialSpectrum(wavelengthSupport, 'LconeIsolating');

[MconeIncrementRadiancePhotons, MconeDecrementRadiancePhotons] = ...
    RGCMosaicConstructor.helper.simulateExperiment.stimulusDifferentialSpectrum(wavelengthSupport, 'MconeIsolating');


% Load the 2-deg Stockman cone fundamentals on wavelength support matching the display
coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), wavelengthSupport);

% Weighting spectra: Stockman-Sharpe 2 deg cone fundamentals
theLconeWeightingSpectrum = squeeze(coneFundamentals(:,1));
theMconeWeightingSpectrum = squeeze(coneFundamentals(:,2));
theLconeWeightingSpectrum = theLconeWeightingSpectrum / sum(theLconeWeightingSpectrum(:));
theMconeWeightingSpectrum = theMconeWeightingSpectrum / sum(theMconeWeightingSpectrum(:));

% Compute L/M-cone fundamental weighted PSFs
nRows = size(thePSFs{1,1}.data,1);
nCols = size(thePSFs{1,1}.data,2);
theLconeWeightedPSFs = cell(numel(sampledXpos), numel(sampledYpos));
theMconeWeightedPSFs = cell(numel(sampledXpos), numel(sampledYpos));

maxWeightedPSFs = 0;
for ix = 1:numel(sampledXpos)
    for iy = 1:numel(sampledYpos)
        thePSFatCurrentEccentricity = thePSFs{ix,iy};
        nRows = size(thePSFatCurrentEccentricity.data,1);
        nCols = size(thePSFatCurrentEccentricity.data,2);
        theCurrentLconeWeightedPSF = zeros(nRows, nCols);
        theCurrentMconeWeightedPSF = zeros(nRows, nCols);
        for iWave = 1:size(thePSFs{1,1}.data,3)
            if (thePSFatCurrentEccentricity.supportWavelength(iWave) ~= wavelengthSupport(iWave))
                error('Wavelengths not matched')
            end
            thePSFslice = squeeze(thePSFatCurrentEccentricity.data(:,:,iWave));
            theCurrentLconeWeightedPSF = theCurrentLconeWeightedPSF + theLconeWeightingSpectrum(iWave) * thePSFslice;
            theCurrentMconeWeightedPSF = theCurrentMconeWeightedPSF + theMconeWeightingSpectrum(iWave) * thePSFslice;
        end % iWave
        theLconeWeightedPSFs{ix,iy} = theCurrentLconeWeightedPSF;
        theMconeWeightedPSFs{ix,iy} = theCurrentMconeWeightedPSF;
        if (max(theCurrentLconeWeightedPSF(:)) > maxWeightedPSFs)
            maxWeightedPSFs = max(theCurrentLconeWeightedPSF(:));
        end
        if (max(theCurrentMconeWeightedPSF(:)) > maxWeightedPSFs)
            maxWeightedPSFs = max(theCurrentMconeWeightedPSF(:));
        end
    end  %iY
end % iX


hFig = figure(2222); clf;
set(hFig, 'Color', [1 1 1]);

subplot(2,1,1)
plot(wavelengthSupport, LconeIncrementRadiancePhotons, 'ro-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.7 0 0.], 'Color', [0.7 0 0.]);
hold on;
plot(wavelengthSupport, MconeIncrementRadiancePhotons, 'go-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.9 0.5], 'MarkerEdgeColor', [0 0.6 0], 'Color', [0 0.6 0]);
plot(wavelengthSupport, achromaticIncrementRadiancePhotons, 'ko-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerEdgeColor', [0.3 0.3 0.3]);
grid on; box off
legend({'L-increment', 'M-increment', 'light-increment'})

set(gca, 'XLim', [400 700], 'FontSize', 16)
xlabel('wavelength (nm)');
ylabel('energy(stim-background)');
title('stimulus spectral power distribution')

subplot(2,1,2)
plot(wavelengthSupport, coneFundamentals(:,1), 'ro-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.7 0 0.], 'Color', [0.7 0 0.]);
hold on;
plot(wavelengthSupport, coneFundamentals(:,2), 'go-', 'LineWidth', 1.5, ...
    'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.9 0.5], 'MarkerEdgeColor', [0 0.6 0], ...
    'Color', [0 0.6 0]);
legend('L-cone', 'M-cone')
set(gca, 'XLim', [400 700], 'FontSize', 16);
grid on; box off
ylabel('sensitivity');
xlabel('wavelength (nm)');
title('cone fundamentals (foveal)');




hFig = figure(2223); clf;
ax = subplot(1,2,1);
imagesc(ax, sampledXpos, sampledYpos, theOptimalStrehlRatio');
colorbar(ax)
title('optimal Strehl ratio')
axis 'image'; axis 'xy';

ax = subplot(1,2,2);
imagesc(ax, sampledXpos, sampledYpos, theOptimalStrehlRatioDefocusDiopters');
colorbar(ax)
title('optimal defocus (D)')
axis 'image'; axis 'xy';



thePSFatCurrentEccentricity = thePSFs{1,1};
wavelengthSupport = thePSFatCurrentEccentricity.supportWavelength;
visualizedWavelengths = 440:10:680;

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', numel(sampledYpos), ...
        'colsNum', numel(sampledXpos), ...
        'heightMargin',  0.01, ...
        'widthMargin',    0.01, ...
        'leftMargin',     0.01, ...
        'rightMargin',    0.0, ...
        'bottomMargin',   0.00, ...
        'topMargin',      0.00);

if (1==2)
% Setup figure
hFig = figure(3); clf;
set(hFig, 'Position', [10 10 2000 1100]);

if (numel(visualizedWavelengths)>1)
    videoFilename = strrep(matFile, '.mat', '');
    videoOBJ = VideoWriter(videoFilename, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();
end

includeStrehlRatiosForEachWavelength = true;
if (includeStrehlRatiosForEachWavelength)
    % AO PSF at 550
    opticsParamsCopy = opticsParams;
    opticsParamsCopy.modification = 'adaptiveOptics6MM';
    opticsParamsCopy.refractiveErrorDiopters =  0.0;

    [~,idx550nm] = min(abs(wavelengthSupport - 550));
    for ix = 1:numel(sampledXpos)
        for iy = 1:numel(sampledYpos)
            eccentricityPosDegs = [sampledXpos(ix) sampledYpos(iy)];
            theConeMosaic = theConeMosaics{ix,iy};
            [~, theAOPSF] = RGCMosaicConstructor.helper.optics.generate(...
                    theConeMosaic, eccentricityPosDegs, opticsParamsCopy, ...
                    'visualizeStrehlRatioOptimization', false);
            theAOPSF = squeeze(theAOPSF.data(:,:,idx550nm));
            theAOPSFpeakAt550(ix,iy) = max(theAOPSF(:));
        end
    end
end


for idx = 1:numel(visualizedWavelengths)
    [~,visualizedWavelengthIndex] = min(abs(wavelengthSupport - visualizedWavelengths(idx)));
    theVisualizedWavelength = wavelengthSupport(visualizedWavelengthIndex);

    for ix = 1:numel(sampledXpos)
        for iy = 1:numel(sampledYpos)
            eccentricityPosDegs = [sampledXpos(ix) sampledYpos(iy)];
            
            thePSFatCurrentEccentricity = thePSFs{ix,iy};
            theConeMosaic = theConeMosaics{ix,iy};


            % Generate the PSF/cone mosaic plot for the current eccentricityPosDegs and visualizedWavelength
            if ((ix == 5) && (iy == 1))
                axPos = subplotPosVectors(iy,ix).v;
                axPos(2) = axPos(2) + axPos(4)*0.13;
                axPos(4) = axPos(4)*0.85;
                ax = subplot('Position', axPos);
                baseline = 0; faceAlpha = 0.4; lineWidth = 1.0;

                cla(ax);
                hold(ax, 'on');
                faceColor = [1 0.2 0.4]; edgeColor = [1 0.2 0.4]*0.5;
                RGCMosaicConstructor.visualize.xyDataAsShadedArea(ax,...
                    wavelengthSupport, coneFundamentals(:,1)/max(coneFundamentals(:)), ...
                    baseline, faceColor, edgeColor, faceAlpha, lineWidth);
                

                faceColor = [0.2 1.0 0.5]; edgeColor = [0.2 1.0 0.5]*0.5;
                RGCMosaicConstructor.visualize.xyDataAsShadedArea(ax,...
                    wavelengthSupport, coneFundamentals(:,2)/max(coneFundamentals(:)), ...
                    baseline, faceColor, edgeColor, faceAlpha, lineWidth);

                plot(ax, wavelengthSupport(visualizedWavelengthIndex)*[1 1], [0 1], '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
                

                faceColor = [1 0.2 0.4]; edgeColor = [1 0.2 0.4]*0.5;
                scatter(ax, wavelengthSupport(visualizedWavelengthIndex), coneFundamentals(visualizedWavelengthIndex,1)/max(coneFundamentals(:)), 14*14, ...
                    'LineWidth', 1.5, 'MarkerFaceAlpha', 0.8, ...
                    'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);

                faceColor = [0.2 1.0 0.5]; edgeColor = [0.2 1.0 0.5]*0.5;
                scatter(ax, wavelengthSupport(visualizedWavelengthIndex), coneFundamentals(visualizedWavelengthIndex,2)/max(coneFundamentals(:)), 14*14, ...
                    'LineWidth', 1.5, 'MarkerFaceAlpha', 0.8, ...
                    'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);

                grid(ax, 'on'); box(ax, 'on');
                %xlabel(ax, 'wavelength (nm)');
                xTicks = [theVisualizedWavelength];
                set(ax, 'FontSize', 20, 'XLim', [wavelengthSupport(1) wavelengthSupport(end)], ...
                    'XLim', [390 700], 'XTick', xTicks, 'XTickLabel', sprintf('%d nm\n',xTicks), ...
                    'YLim', [0 1], 'YTick', [], 'YTickLabel', {});
                xtickangle(ax, 0);
                hold(ax, 'off');
            else
                ax = subplot('Position', subplotPosVectors(iy,ix).v);
                theTitle = sprintf('%d nm (Sr: %2.3f)', ...
                        thePSFatCurrentEccentricity.supportWavelength(visualizedWavelengthIndex), theOptimalStrehlRatio(ix,iy));
                theTitle = '';

                if (includeStrehlRatiosForEachWavelength)
                    theStrehlRatioAtThisWavelength = ...
                        max(max(squeeze(thePSFatCurrentEccentricity.data(:,:,visualizedWavelengthIndex)))) / theAOPSFpeakAt550(ix,iy);
                    infoString = sprintf('%.2f', theStrehlRatioAtThisWavelength);
                else
                    infoString = '';
                end

                RGCMosaicConstructor.helper.optics.visualizePSFAtWavelengthOnTopOfConeMosaic(...
                    hFig, ax, theConeMosaic, thePSFatCurrentEccentricity, ...
                    visualizedWavelengthIndex, max(thePSFatCurrentEccentricity.data(:)), ...
                    theTitle, ...
                    'XLimsArcMin', [-3 3], ...
                    'YLimsArcMin', [-3 3], ...
                    'XTicksArcMin', -4:0:4, ...
                    'YTicksArcMin', -4:0:4, ...
                    'withInfoString', infoString);
            end      
        end
    end

    if (numel(visualizedWavelengths)>1)
        % Write frame data
        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end
end % for inFocusWavelengthIndex


if (numel(visualizedWavelengths)>1)
    % Save video
    videoOBJ.close()
end

end

% Generate figure with L-cone weighted PSFs as a function of eccentricity
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 2000 1100], 'Name', 'LconePSFs');
for ix = 1:numel(sampledXpos)
    for iy = 1:numel(sampledYpos)

        theConeMosaic = theConeMosaics{ix,iy};
        theConePSFatCurrentEccentricity = thePSFs{ix,iy};
        theConeWeightedPSF = theLconeWeightedPSFs{ix,iy};
        nRows = size(theConeWeightedPSF,1);
        nCols = size(theConeWeightedPSF,2);
         % Replace the 3D PSF with a 2D cone-weighted PSF
        theConePSFatCurrentEccentricity.data = reshape(theConeWeightedPSF, [nRows nCols 1]);

        % Generate the PSF/cone mosaic plot for the current eccentricityPosDegs and visualizedWavelength
        ax = subplot('Position', subplotPosVectors(iy,ix).v);
        theTitle = sprintf('(x,y)=(%d,%d)', sampledXpos(ix), sampledYpos(iy));
        theTitle = '';
        RGCMosaicConstructor.helper.optics.visualizePSFAtWavelengthOnTopOfConeMosaic(...
            hFig, ax, theConeMosaic, theConePSFatCurrentEccentricity, ...
            1, maxWeightedPSFs, ...
            theTitle, ...
            'XLimsArcMin', [-3 3], 'YLimsArcMin', [-3 3], 'XTicksArcMin', -4:0:4,'YTicksArcMin', -4:0:4);        
    end
end
thePDFfileName = strrep(matFile, '.mat', '_LconePSFs.pdf');
NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);

% Generate figure with M-cone weighted PSFs as a function of eccentricity
hFig = figure(2); clf;
set(hFig, 'Position', [10 10 2000 1100], 'Name', 'MconePSFs');
for ix = 1:numel(sampledXpos)
    for iy = 1:numel(sampledYpos)

        theConeMosaic = theConeMosaics{ix,iy};
        theConePSFatCurrentEccentricity = thePSFs{ix,iy};
        theConeWeightedPSF = theMconeWeightedPSFs{ix,iy};
        nRows = size(theConeWeightedPSF,1);
        nCols = size(theConeWeightedPSF,2);
        % Replace the 3D PSF with a 2D cone-weighted PSF
        theConePSFatCurrentEccentricity.data = reshape(theConeWeightedPSF, [nRows nCols 1]);

        % Generate the PSF/cone mosaic plot for the current eccentricityPosDegs and visualizedWavelength
        ax = subplot('Position', subplotPosVectors(iy,ix).v);
        theTitle = sprintf('(x,y)=(%d,%d)', sampledXpos(ix), sampledYpos(iy));
        theTitle = '';
        RGCMosaicConstructor.helper.optics.visualizePSFAtWavelengthOnTopOfConeMosaic(...
            hFig, ax, theConeMosaic, theConePSFatCurrentEccentricity, ...
            1, maxWeightedPSFs, ...
            theTitle, ...
            'XLimsArcMin', [-3 3], 'YLimsArcMin', [-3 3], 'XTicksArcMin', -4:0:4,'YTicksArcMin', -4:0:4);        
    end
end
thePDFfileName = strrep(matFile, '.mat', '_MconePSFs.pdf');
NicePlot.exportFigToPDF(thePDFfileName, hFig, 300);
