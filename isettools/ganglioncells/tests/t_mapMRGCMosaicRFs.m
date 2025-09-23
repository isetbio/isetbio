% Script to map the RF centers using m-sequence stimuli
%
% Usage:
%{
    t_mapMRGCMosaicRFs     
%}

% Initialize session
close all; clear all;

% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
ASPPManager = AppleSiliconParPoolManager('conservative');
%ASPPManager = AppleSiliconParPoolManager(12);

employRFCenterOverlappingMosaic = true;

% Examined chroma-spatial variance tradeoff values [0: minimize chromatic variance 1: minimize spatial variance]
chromaticSpatialVarianceTradeoff = 1.0;

% Which mRGC/cone mosaic lattice to use
sourceLatticeSizeDegs = 64; whichEye = 'right eye';

% Generate a 18x18 deg RGC mosaic centered at (-12,0) in the right eye
mosaicEccDegs = [-12 0]; mosaicSizeDegs = 17*[1 1];
mosaicEccDegs = [-20 0]; mosaicSizeDegs = 27*[1 1];

% Extra size  
extraSizeDegs = 1;
mosaicSizeDegs = mosaicSizeDegs+extraSizeDegs*[1 1];

% Grid of (X,Y)-positions, (W,H)-sizes on which the surround will be optimized 
% The size vector controls the area over which we will searh for unique # of center conesNum 
rfMappingPatchGrids = [];
%rfMappingPatchGrids(size(rfMappingPatchGrids,1)+1,:) = [-12 0 0.5 0.5];


xo = -7:-2:-31;
xo = -31;
for rfMappingPatchIndex = 1:numel(xo)
    sizeDegs = 0.5 + abs(xo(rfMappingPatchIndex))*0.02;
    rfMappingPatchGrids(size(rfMappingPatchGrids,1)+1,:)  = [xo(rfMappingPatchIndex) 0 sizeDegs sizeDegs];
end


% Use mosaic with RFcenter overlap
operateOnCenterConnectedMosaicsWithOverlap = true;

% Actions to perform
% Input cone mosaic responses
performComputeInputConeMosaicResponsesAction = ~true;
performComputeMRGCMosaicResponsesAction = ~true;
performVisualizeMRGCMosaicRFmaps = true;

% Flag indicating whether to generate a summary figure over all positions examined
generateSummaryRFmapsOverAllPositionsFigure = true;
includeSummaryDataFromEJpapers = true;
summaryDataEccUnits = 'MMs';  % choose between {'degs', 'MMs'}

% Flag indicating whether to visualize the input cone mosaic responses
visualizeSampleInputConeMosaicResponses = ~true;

employAdaptiveOptics = true;

% Normal subject with optimized Strehl ratio PSF
opticsParamsStruct = struct(...
    'ZernikeDataBase', 'Polans2015', ...
    'subjectID', 6, ...
    'pupilDiameterMM', 3.0, ...
    'noLCA', false, ...
    'zeroCenterPSF', true, ...
    'modification', 'optimizedStrehlRatio' ... % choose from {'none', 'adaptiveOptics', subtractedCentralRefraction', 'optimizedStrehlRatio', etc}
);

if (employAdaptiveOptics)
    % Adaptive optics with 6mm pupil
    opticsParamsStruct = struct(...
        'ZernikeDataBase', 'Polans2015', ...
        'subjectID', 6, ...
        'pupilDiameterMM', 6.0, ...
        'noLCA', true, ...
        'zeroCenterPSF', true, ...
        'modification', 'adaptiveOptics' ... % choose from {'none', 'adaptiveOptics', subtractedCentralRefraction', 'optimizedStrehlRatio', etc}
    );
end

RFmappingParamsStruct = struct(...
    'backgroundChromaticity', [0.301 0.301], ...
    'backgroundLuminanceCdM2', 50.0, ...
    'chromaticity', 'Achromatic', ...
    'coneFundamentalsOptimizedForStimPosition', false, ...
    'ternaryInsteadOfBinaryMsequence', false, ...
    'rfPixelsAcross', 75, ...
    'mSequenceBitLength', 13, ...
    'frameBatchSize', 256, ...
    'resolutionDegs', [], ...                       % to be determined separately for each optimization position
    'positionDegs',[], ... 
    'sizeDegs', [] ... 
);

mRGCMosaicParamsStruct = struct(...
    'whichEye', whichEye, ...
    'eccentricityDegs',  mosaicEccDegs, ...
    'sizeDegs', mosaicSizeDegs, ...
    'chromaticSpatialVarianceTradeoff', chromaticSpatialVarianceTradeoff, ...
    'rfCentersOverlap', operateOnCenterConnectedMosaicsWithOverlap);


% Run all patches
for rfMappingPatchIndex = 1:size(rfMappingPatchGrids,1)
    % Run this position
    RFmappingParamsStruct.positionDegs =  rfMappingPatchGrids(rfMappingPatchIndex,1:2);
    RFmappingParamsStruct.sizeDegs =  rfMappingPatchGrids(rfMappingPatchIndex,3:4);
    
    [theMinorSigmasDegs, theMajorSigmasDegs, theTemporalEquivalentEccentricitiesDegs, ...
     theMinorSigmasMicrons, theMajorSigmasMicrons, theTemporalEquivalentEccentricitiesMMs] = ...
        RGCMosaicConstructor.helper.simulateExperiment.mSequenceRFmapping(...
            mRGCMosaicParamsStruct, opticsParamsStruct, RFmappingParamsStruct, ...
            performComputeInputConeMosaicResponsesAction, ...
            performComputeMRGCMosaicResponsesAction, ...
            performVisualizeMRGCMosaicRFmaps, ...
            'visualizeInputConeMosaicResponses', visualizeSampleInputConeMosaicResponses);

    if (generateSummaryRFmapsOverAllPositionsFigure)
        micronsPerDegreeInMacaqueRetina = 222;

        % Initialize figure
        if (rfMappingPatchIndex == 1)
            figNo = 500;
            hFig = figure(figNo); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure')
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)]);
            ax = theAxes{1,1};
        end
    
        % Add data from current mosaic
        % Compute radial eccentricities
        theTemporalEquivalentRadiiDegs = sqrt(sum(theTemporalEquivalentEccentricitiesDegs.^2,2));
        theTemporalEquivalentRadiiMMs = sqrt(sum(theTemporalEquivalentEccentricitiesMMs.^2,2));

        if (strcmp(summaryDataEccUnits,'MMs'))
            ecc = theTemporalEquivalentRadiiMMs(:);
        else
            ecc = theTemporalEquivalentRadiiDegs(:);
        end

        hold(ax, 'on');

        % We are plotting RF diameters, so multiply sigma by 2
        scatter(ax, ecc, 2*theMinorSigmasMicrons(:), (ff.markerSize-6)^2,'o', ...
        'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5, ...
        'LineWidth', 0.5);
        scatter(ax, ecc, 2*theMajorSigmasMicrons(:), (ff.markerSize-6)^2,'o', ...
            'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1], 'MarkerFaceAlpha', 0.5, ...
            'LineWidth', 0.5);
        ylabel(ax, 'RF diameter (2 \times \sigma), microns');

        set(ax, 'YLim', [0 80], 'YTick', 0:10:100);
        if (strcmp(summaryDataEccUnits,'MMs'))
            set(ax, 'XLim', [1 9], 'XTick', 1:10);
            xlabel(ax, 'temporal equivalent eccentricity, mm')
        else
            XLimsDegs = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', [1 9]);
            set(ax, 'XLim', [floor(XLimsDegs(1)) ceil(XLimsDegs(2))], 'XTick', 0:5:100);
            xlabel(ax, 'temporal equivalent eccentricity, degs')
        end
    

        % Include EJ's data
        if ((includeSummaryDataFromEJpapers) && (rfMappingPatchIndex == size(rfMappingPatchGrids,1)))
            dataOut = RGCMosaicConstructor.publicData.ChichilniskyLab.rfDiameters();
            d = dataOut('3.5mm');
            eccMM = abs(d.patchEccMM) + -0.05+d.randomJitter*0.1;
            eccDegs = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', eccMM);
            if (strcmp(summaryDataEccUnits,'MMs'))
                ecc = eccMM;
            else
                ecc = eccDegs;
            end

            plot(ax, ecc, d.rfDiameterMicrons, 'g.', 'MarkerSize', 20);
            d = dataOut('8.5mm');
            eccMM = abs(d.patchEccMM) + -0.05+d.randomJitter*0.1;
            eccDegs = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', eccMM);
            if (strcmp(summaryDataEccUnits,'MMs'))
                ecc = eccMM;
            else
                ecc = eccDegs;
            end
            plot(ax, ecc, d.rfDiameterMicrons, 'g.', 'MarkerSize', 20);
            fprintf('Need to load scanned data here\n');
        end


        % Finalize figure
        if (rfMappingPatchIndex  == size(rfMappingPatchGrids,1))
            % Finalize figure using the Publication-Ready format
            PublicationReadyPlotLib.applyFormat(ax,ff);
            PublicationReadyPlotLib.offsetAxes(ax, ff, [], []);

            % Export figure
            theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            pdfExportSubDir = 'summaryPlots';
            if (strcmp(summaryDataEccUnits,'degs'))
                thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'RFdiametersMicrons_vs_eccentricityDegs.pdf');
            else
                thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'RFdiametersMicrons_vs_eccentricityMMs.pdf');
            end
            NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
        end
    end % if (generateSummaryRFmapsOverAllPositionsFigure)
end % for rfMappingPatchIndex


