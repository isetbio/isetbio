function measureMRGCRFcentersAcrossXYeccentricityForDifferentSubjects()

    % Vertical eccentricities examined
    eccX = [-20 -16 -12 -10 -8:1:8 10 12 20];
    eccY = [-8:1:8];

    % Remaining computations
    eccX = [12 20];
    eccY = [-8:1:8];


    % Horizontal eccentricities examined
    eccX = [-7];
    eccY = [0];

    eccX = [-10 -8:1:8 10];
    eccY = [-6:1:6];

    eccX = [0];
    eccY = [-8:1:-1];


    % Optics
    %newZernikeDataBase = 'Artal2012';
    %newSubjectRankOrder = 3;

    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 6;

    newZernikeDataBase = 'Polans2015';
    newSubjectRankOrder = 6;
    newPupilDiamMM = 3.5;
    newWavefrontSpatialSamplesNum = 701;


    maxVisualizedRFs = 64;

    % Generate the midgetRGCmosaic and optics
    generateTheComponents = ~true;
    if (generateTheComponents)
        generateComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY);
        visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs);
    end

    % Modify the optics
    changeOpticsSubject = ~true;
    if (changeOpticsSubject)
        modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, ...
            newSubjectRankOrder, newPupilDiamMM, newWavefrontSpatialSamplesNum, eccX, eccY);
        visualizeComponents(newZernikeDataBase, newSubjectRankOrder, ...
            eccX, eccY, maxVisualizedRFs);
    end
    
    
    % Compute RF maps (center)
    mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesData';

    computeTheVisualRFmaps = true;
    if (computeTheVisualRFmaps)
        computeVisualRF(newZernikeDataBase, newSubjectRankOrder, ...
            eccX, eccY, maxVisualizedRFs, mappedRFsDir);
    end

    % Fit the retinal and visual RF maps
    fitAndVisualizeTheRetinalAndVisualRFmaps = true;
    if (fitAndVisualizeTheRetinalAndVisualRFmaps)
        visualizeFits = ~true;
        fitTheRetinalAndVisualRFmaps(newZernikeDataBase, newSubjectRankOrder, ...
            eccX, eccY, mappedRFsDir, visualizeFits);
    end

    % Summarize analyses
    summarizeAnalyses = ~true;
    if (summarizeAnalyses)
        summarizeAnalyzedData(newZernikeDataBase, newSubjectRankOrder, ...
            eccX, eccY, mappedRFsDir);
    end
end


function summarizeAnalyzedData(ZernikeDataBase, subjectRankOrder, eccX, eccY, mappedRFsDir)

    visualizedFOVdegs = 0.15;

%     hFig = figure(2); clf;
%     set(hFig, 'Position', [10 10 1650 1050], 'Color', [1 1 1]);
%     maxVisualizedRFs = 7;
%     plotTheMidgetRGC(visualizedFOVdegs, ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs, mappedRFsDir, hFig);
%     NicePlot.exportFigToPDF('RetinalRFs.pdf', hFig, 300);

    % The grid of PSFs
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1650 1050], 'Color', [1 1 1]);
    plotThePSFgrid(visualizedFOVdegs, ZernikeDataBase, subjectRankOrder, eccX, eccY, mappedRFsDir);
    NicePlot.exportFigToPDF('PSFs.pdf', hFig, 300);

    

    plotTheRcGrid(mappedRFsDir);

end


function plotTheMidgetRGC(visualizedFOVdegs, ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs, mappedRFsDir, hFig)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    mosaicEccDegsGrid = [eccXGrid eccYGrid];

    rowsNum = numel(eccY);
    colsNum = numel(eccX);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', rowsNum, ...
       'colsNum', colsNum, ...
       'heightMargin',  0.005, ...
       'widthMargin',    0.005, ...
       'leftMargin',     0.011, ...
       'rightMargin',    0.005, ...
       'bottomMargin',   0.005, ...
       'topMargin',      0.005);

    coVisualizePSF = false;
    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = mosaicEccDegsGrid(iPos,:);
 
        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        load(fName, 'thePSFData', 'theMidgetRGCmosaic');

        subplotCol = floor((iPos-1)/rowsNum)+1;
        subplotRow = rowsNum - mod(iPos-1,rowsNum);

        axesHandle = subplot('Position', subplotPosVectors(subplotRow,subplotCol).v);

        if (coVisualizePSF)
            if (isfield(thePSFData, 'vLambdaWeightedData'))
                psfData.supportXdegs = thePSFData.psfSupportXdegs;
                psfData.supportYdegs = thePSFData.psfSupportYdegs;
                psfData.data = thePSFData.vLambdaWeightedData;
            else
                psfData.supportXdegs = thePSFData.psfSupportXdegs;
                psfData.supportYdegs = thePSFData.psfSupportYdegs;
                psfData.data = thePSFData.unrotatedData/max(thePSFData.unrotatedData(:));
            end
    
            theMidgetRGCmosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', axesHandle, ...
                'maxVisualizedRFs', maxVisualizedRFs, ...
                'retinalMeridianAxesLabeling', false, ...
                'withSuperimposedPSF', psfData, ...
                'xRangeDegs', visualizedFOVdegs, ...
                'yRangeDegs', visualizedFOVdegs, ...
                'noXLabel', true, ...
                'noYLabel', true, ...
                'fontSize', 12);
        else
            theMidgetRGCmosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', axesHandle, ...
                'maxVisualizedRFs', maxVisualizedRFs, ...
                'retinalMeridianAxesLabeling', false, ...
                'xRangeDegs', visualizedFOVdegs, ...
                'yRangeDegs', visualizedFOVdegs, ...
                'noXLabel', true, ...
                'noYLabel', true, ...
                'fontSize', 12);
        end

        set(axesHandle,'XTickLabel', {}, 'YTickLabel', {});

        if (subplotRow == 1)
            title(axesHandle, sprintf('x:%2.0f degs', mosaicEccDegs(1)), 'FontWeight', 'normal');
        end
        if (subplotCol == 1)
            ylabel(axesHandle, sprintf('y:%2.0f degs', mosaicEccDegs(2)));
        end
        drawnow;
    end

end


function plotThePSFgrid(visualizedFOVdegs, ZernikeDataBase, subjectRankOrder, eccX, eccY, mappedRFsDir)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    mosaicEccDegsGrid = [eccXGrid eccYGrid];

    rowsNum = numel(eccY);
    colsNum = numel(eccX);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', rowsNum, ...
       'colsNum', colsNum, ...
       'heightMargin',  0.005, ...
       'widthMargin',    0.005, ...
       'leftMargin',     0.011, ...
       'rightMargin',    0.005, ...
       'bottomMargin',   0.005, ...
       'topMargin',      0.005);
    

    % Params for PSF contour plot
    % Color map
    cmap = brewermap(1024,'blues');
    
    % Transparency level
    alpha = 0.75;
    
    % Color of the contour lines
    contourLineColor = [0.2 0.2 0.2];
    
    % Levels of the contour lines
    zLevels = 0.1:0.1:0.95;

    
    for iPos = 1:numel(eccXGrid)
        % The position analyzed
        mosaicEccDegs = mosaicEccDegsGrid(iPos,:);
 
        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        load(fName, 'thePSFData', 'theMidgetRGCmosaic');

        % Load the computed RF maps data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        fNameRFmaps = fullfile(mappedRFsDir,fNameRFmaps);
        load(fNameRFmaps, 'visualRFspatialSupportDegs', ...
                          'retinalRFcenterMaps', 'theVisualRFmaps', 'rgcIndicesOfAnalyzedRFs');


        subplotCol = floor((iPos-1)/rowsNum)+1;
        subplotRow = rowsNum - mod(iPos-1,rowsNum);
        axesHandle = subplot('Position', subplotPosVectors(subplotRow,subplotCol).v);

        if (isfield(thePSFData, 'vLambdaWeightedData'))
            cMosaic.semiTransparentContourPlot(axesHandle, ...
                 thePSFData.psfSupportXdegs,...
                 thePSFData.psfSupportXdegs, ...
                 thePSFData.vLambdaWeightedData/max(thePSFData.vLambdaWeightedData(:)), ...
                 zLevels, cmap, alpha, contourLineColor, ...
                'lineWidth', 1.0);
        else
            cMosaic.semiTransparentContourPlot(axesHandle, ...
                 thePSFData.psfSupportXdegs,...
                 thePSFData.psfSupportXdegs, ...
                 thePSFData.unrotatedData/max(thePSFData.unrotatedData(:)), ...
                 zLevels, cmap, alpha, contourLineColor, ...
                'lineWidth', 1.0);
        end

        set(axesHandle, 'xLim', 0.5*visualizedFOVdegs*[-1 1], ...
                        'yLim', 0.5*visualizedFOVdegs*[-1 1], ...
                        'xTick', -0.15:0.025:0.15, ...
                        'yTick', -0.15:0.025:0.15, ...
                        'XTickLabel', {}, 'YTickLabel', {}, ...
                        'fontSize', 12);
        axis(axesHandle, 'square');
        axis(axesHandle, 'xy');
        grid(axesHandle, 'on');
        box(axesHandle, 'on');
       
        if (subplotRow == 1)
            title(axesHandle, sprintf('x:%2.0f degs', mosaicEccDegs(1)), 'FontWeight', 'normal');
        end
        if (subplotCol == 1)
            ylabel(axesHandle, sprintf('y:%2.0f degs', mosaicEccDegs(2)));
        end
        drawnow;

    end
end

function plotTheRcGrid(mappedRFsDir)

     % Load the summary data
     load(fullfile(mappedRFsDir, 'summaryData.mat'), 'mosaicEccDegsGrid', 'visualRc', 'retinalRc');

     hFig = figure(3);
     clf;
     set(hFig, 'Position', [10 10 800 800], 'Color', [1 1 1]);


     % Determine the rance of Rcs
     minRc = Inf; maxRc = 0;
     for iPos = 1:size(mosaicEccDegsGrid,1)
         mosaicEccDegs = mosaicEccDegsGrid(iPos,:);
         theVisualRcs = visualRc{iPos};
         theRetinalRcs = retinalRc{iPos};

         if (~isempty(theVisualRcs))
             m = mean(theVisualRcs);
             if (m< minRc)
                 minRc = m;
             end
             if (m> maxRc)
                 maxRc = m;
             end
             m = mean(theRetinalRcs);
             if (m< minRc)
                 minRc = m;
             end
             if (m> maxRc)
                 maxRc = m;
             end
         end
     end


     minMarkerSize = 5;
     maxMarkerSize = 20;

     for iPos = 1:size(mosaicEccDegsGrid,1)
         mosaicEccDegs = mosaicEccDegsGrid(iPos,:);
         theVisualRcs = visualRc{iPos};
         theRetinalRcs = retinalRc{iPos};

         if (~isempty(theVisualRcs))
             meanVisualRc = mean(theVisualRcs);
             meanRetinalRc = mean(theRetinalRcs);

             % the visual Rc
             markerSize = minMarkerSize^2 +(maxMarkerSize)^2 * (meanVisualRc-minRc)/(maxRc-minRc);
             scatter(mosaicEccDegs(1), mosaicEccDegs(2), markerSize, ...
                 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0], ...
                 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
                 
             hold on; 

             % the retinal Rc
             markerSize = minMarkerSize^2 + (maxMarkerSize)^2 * (meanRetinalRc-minRc)/(maxRc-minRc);
             scatter(mosaicEccDegs(1), mosaicEccDegs(2), markerSize, ...
                 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1], ...
                 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

         end

     end

     markerSize = minMarkerSize^2;
     scatter(-3.5, -9.5, markerSize, ...
                 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0 0 0], ...
                 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
     text(-3.0, -9.5, sprintf('min Rc: %2.1f arc min', minRc), 'FontSize', 12);

     markerSize = minMarkerSize^2 + (maxMarkerSize)^2;
     scatter(3.5, -9.5, markerSize, ...
                 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0 0 0], ...
                 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
     text(4.0, -9.5, sprintf('max Rc: %2.1f arc min', maxRc), 'FontSize', 12);

     axis 'square';
     grid on;
     box off;
     set(gca, 'XLim', [-10 10], 'YLim', [-10 10], 'XTick', -10:1:10, 'YTick', -10:1:10);
     set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0);
     set(gca, 'TickDir', 'both');
     set(gca, 'FontSize', 16);
     xlabel('space (degs)');
     ylabel('space (degs)');



end

function fitTheRetinalAndVisualRFmaps(ZernikeDataBase, subjectRankOrder, eccX, eccY, mappedRFsDir, visualizeFits)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    mosaicEccDegsGrid = [eccXGrid eccYGrid];
    retinalRc = cell(1, numel(eccXGrid));
    visualRc = cell(1, numel(eccXGrid));

    for iPos = 1:numel(eccXGrid)

        fprintf('\nanalyzing position %d of %d\n', iPos, numel(eccXGrid));
        % The position analyzed
        mosaicEccDegs = mosaicEccDegsGrid(iPos,:);

        % Load data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        fNameRFmaps = fullfile(mappedRFsDir,fNameRFmaps);

        load(fNameRFmaps, 'theMidgetRGCmosaic', 'opticsParams', 'visualRFspatialSupportDegs', ...
                          'retinalRFcenterMaps', 'theVisualRFmaps', 'rgcIndicesOfAnalyzedRFs');

        if (visualizeFits)
            [~, psfEnsemble] = theMidgetRGCmosaic.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
            thePSFData = psfEnsemble{1};
            idx = find(thePSFData.supportWavelength == 550);
            thePSF = thePSFData.data(:,:,idx);
            
            figure(5); clf;
            ax = subplot(2,3,1);
            imagesc(ax, thePSFData.supportX/60, thePSFData.supportY/60, thePSF/max(thePSF(:)));
            axis(ax,'image'); axis(ax, 'xy');
            set(ax, 'XLim', 0.25*[-1 1], 'YLim', 0.25*[-1 1]);
            title(ax, 'PSF');
    
            xLims = theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1) + [-0.25 0.25];
            yLims = theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2) + [-0.25 0.25];
        end

        % Initialize memory
        thePositionRetinalRc = zeros(1, numel(rgcIndicesOfAnalyzedRFs));
        thePositionVisualRc = zeros(1, numel(rgcIndicesOfAnalyzedRFs));

        parfor iRGC = 1:numel(rgcIndicesOfAnalyzedRFs)
            % Fit the retinal RFmap
            r = retinalRFcenterMaps{iRGC};
            centroidDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs(r.inputConeIndices,:),1);
            theRetinalRFmap = r.centerRF;
            theFittedRetinalRFmap = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                r.spatialSupportDegsX, r.spatialSupportDegsY, ...
                theRetinalRFmap, ...
                'flatTopGaussian', true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/2 2/1], ...
                'forcedCentroidXYpos', centroidDegs , ...
                'globalSearch', true, ...
                'multiStartsNum', 8, ...
                'useParallel', false);
            thePositionRetinalRc(iRGC) = sqrt(prod(theFittedRetinalRFmap.characteristicRadii))*60;

            % Fit the visual RFMap
            theVisualRFmap = double(squeeze(theVisualRFmaps(iRGC,:,:)));
            theVisualRFmap = theVisualRFmap / max(theVisualRFmap(:));

            theFittedVisualRFmap = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                visualRFspatialSupportDegs, visualRFspatialSupportDegs, ...
                theVisualRFmap, ...
                'flatTopGaussian', false, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/2 2/1], ...
                'globalSearch', true, ...
                'multiStartsNum', 8, ...
                'useParallel', false);
            thePositionVisualRc(iRGC) = sqrt(prod(theFittedVisualRFmap.characteristicRadii))*60;

            if (visualizeFits)
                % The retinal RFmap
                ax = subplot(2,3,2);
                imagesc(ax,r.spatialSupportDegsX, r.spatialSupportDegsY, theRetinalRFmap);
                axis(ax,'image'); axis(ax, 'xy');
                set(ax, 'XLim', xLims, 'YLim', yLims);
                title(ax, 'Retinal RF');

                % The fitted retinal RFmap
                ax = subplot(2,3,3);
                imagesc(ax,r.spatialSupportDegsX, r.spatialSupportDegsY, theFittedRetinalRFmap.ellipsoidMap);
                axis(ax,'image'); axis(ax, 'xy');
                set(ax, 'XLim', xLims, 'YLim', yLims);
                title(ax, sprintf('Fitted Retinal RF\n(Rc: %2.3f arc min)', thePositionRetinalRc(iRGC)));

                % The measured and the fitted visual RFMap
                max1 = max(theVisualRFmap(:));
                max2 = max(theFittedVisualRFmap.ellipsoidMap(:));
                maxVisualRF = max([max1 max2]);

                % Measured
                ax = subplot(2,3,4);
                imagesc(ax,visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                           visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                           theVisualRFmap);
                axis(ax,'image'); axis(ax, 'xy');
                set(ax, 'XLim', xLims, 'YLim', yLims, 'Color', [0 0 0], 'CLim', [0 maxVisualRF]);
                title(ax, sprintf('visualRF %s_rank%d', ZernikeDataBase, subjectRankOrder));

                % Fitted
                ax = subplot(2,3,5);
                imagesc(ax,visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                           visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                           theFittedVisualRFmap.ellipsoidMap);
                axis(ax,'image'); axis(ax, 'xy');
                set(ax, 'XLim', xLims, 'YLim', yLims, 'Color', [0 0 0], 'CLim', [0 maxVisualRF]);
                title(ax, sprintf('Fitted VisualRF\n(Rc: %2.3f arc min)', thePositionVisualRc(iRGC)));

                % Residual
                ax = subplot(2,3,6);
                imagesc(ax,visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                           visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                           theVisualRFmap-theFittedVisualRFmap.ellipsoidMap );
                axis(ax,'image'); axis(ax, 'xy');
                set(ax, 'XLim', xLims, 'YLim', yLims, 'Color', [0 0 0], 'CLim', maxVisualRF*[-1 1]);
                title(ax, sprintf('residual visual RF'));
                drawnow
            end % visualizeFits

        end % iRGC

         % Save the data
        save(fullfile(mappedRFsDir, 'summaryData.mat'), 'mosaicEccDegsGrid', 'visualRc', 'retinalRc')

        % Save the visual and retinal Rcs for these RGCs
        retinalRc{iPos} = thePositionRetinalRc;
        visualRc{iPos} = thePositionVisualRc;
    end % iPos

end


function computeVisualRF(ZernikeDataBase, subjectRankOrder, eccX, eccY, centerMostRGCsNumToAnalyze, mappedRFsDir)
    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    opticalImagePositionDegs = 'mosaic-centered';  % [0.05 -0.03];  %'mosaic-centered'

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'theMidgetRGCmosaic', 'opticsParams');

        % Compute the selected subject optics using the saved opticsParams
        % Change pupil diameter to 3.5 mm (to match Johannes experiment)
        opticsParams.pupilDiameterMM = 3.0;

        oiEnsemble = theMidgetRGCmosaic.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                    'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                    'subjectID', opticsParams.testSubjectID, ...
                    'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                    'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                    'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                    'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                    'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);
        theSubjectOptics = oiEnsemble{1};

        % Compute visual RF maps using subspace rev corr
        % Generate a presentation display with a desired resolution
        stimSizeDegs = 0.5*max(theMidgetRGCmosaic.sizeDegs);
        pixelsNum = 256;
        retinalImageResolutionDegs = stimSizeDegs/pixelsNum;
        viewingDistanceMeters = 4;
        theDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            theMidgetRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

        % Stim params for the RF mapping
        stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', [1 1 1], ...
            'contrast', 0.75, ...
            'pixelSizeDegs', retinalImageResolutionDegs, ...
            'stimSizeDegs', stimSizeDegs, ...  % make the stimulus size = 1/2 * RGC mosaic FoV
            'wavelengthSupport', displayGet(theDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(theDisplay, 'viewing distance') ...
            );

        % Hartley (RF mapping) spatial patterns
        fprintf('Generating Hartley patterns\n');
        omega = 11;
        % Compute spatial modulation patterns for the Hartley set
        HartleySpatialModulationPatterns = ...
            rfMappingStimulusGenerator.HartleyModulationPatterns(...
            omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    

        
        % Find the indices of the centerMostRGCsNumToAnalyze RGCs
        relativeRGCpositions = bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs);
        radii = sum(relativeRGCpositions.^2,2);
        [~,sortedRGCindices] = sort(radii, 'ascend');
        rgcsNum = size(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
        rgcIndicesOfAnalyzedRFs = sortedRGCindices(1:min([centerMostRGCsNumToAnalyze, rgcsNum]));

        % Preallocate memory
        nStim = size(HartleySpatialModulationPatterns,1);
        pixelsNum = size(HartleySpatialModulationPatterns,2);
        theVisualRFmaps = zeros(numel(rgcIndicesOfAnalyzedRFs), pixelsNum, pixelsNum, 'single');

        % Compute theNullStimulusScene and the visual RF spatial support
        [~, theNullStimulusScene, visualRFspatialSupportDegs] = ...
            rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, HartleySpatialModulationPatterns(1,:,:), ...
                    'validateScenes', false);


        
        % Visualize first frame
%         if (~isempty(visualizedOpticalImageNum))
%                 
%                 hFig = figure(75);
% 
%                 % Visualize the input cone mosaic with the OI
%                 ax = subplot(2,2,1);
%                 theMidgetRGCmosaic.inputConeMosaic.visualize(...
%                     'figureHandle', hFig, ...
%                     'axesHandle', ax, ...
%                     'withSuperimposedOpticalImage', theOpticalImage);
% 
%                 % Visualize the input cone mosaic response to the OI
%                 ax = subplot(2,2,2);
%                 theMidgetRGCmosaic.inputConeMosaic.visualize(...
%                     'figureHandle', hFig, ...
%                     'axesHandle', ax, ...
%                     'activation', theNoiseFreeConeAbsorptionsCount , ...
%                     'labelCones', false);
% 
%                 % Visualize the midgetRCCmosaic with the OI
%                 ax = subplot(2,2,3);
%                 theMidgetRGCmosaic.visualize(...
%                     'figureHandle', hFig, ...
%                     'axesHandle', ax, ...
%                     'xRangeDegs', theMidgetRGCmosaic.sizeDegs(1), ...
%                     'yRangeDegs', theMidgetRGCmosaic.sizeDegs(2), ...
%                     'withSuperimposedOpticalImage', theOpticalImage);
%                 drawnow;
%         end

        diffResponse = zeros(numel(rgcIndicesOfAnalyzedRFs), nStim);

        parfor iFrame = 1:nStim
            fprintf('Computing responses and RFs for stimulus %d/%d\n', iFrame, nStim);

            theHartleyPattern = HartleySpatialModulationPatterns(iFrame,:,:);

            % Generate scene for the forward Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, theHartleyPattern, ...
                    'validateScenes', false);

            % Compute the mosaic's response to the positive polarity frame
            rPlus = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'opticalImagePositionDegs', opticalImagePositionDegs, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Generate scene for the inverse Hartley pattern
            theRFMappingStimulusFrameScene = ...
                rfMappingStimulusGenerator.generateStimulusMappingFramesOnPresentationDisplay(...
                    theDisplay, stimParams, -theHartleyPattern, ...
                    'validateScenes', false);

            % Compute the mosaic's response to the reverse polarity frames
            rMinus = theMidgetRGCmosaic.compute(...
                    theRFMappingStimulusFrameScene{1}, ...
                    'withOptics', theSubjectOptics, ...
                    'nTrials', 1, ...
                    'opticalImagePositionDegs', opticalImagePositionDegs, ...
                    'theNullScene', theNullStimulusScene, ...
                    'normalizeConeResponsesWithRespectToNullScene', true);

            % Save differentual response to this stimulus frame
            dR = rPlus-rMinus;
            diffResponse(:,iFrame) = squeeze(dR(rgcIndicesOfAnalyzedRFs));
        end

        

       for iFrame = 1:nStim

            theHartleyPattern = HartleySpatialModulationPatterns(iFrame,:,:);
            % Accumulate the RFs
            parfor iRGCindex = 1:numel(rgcIndicesOfAnalyzedRFs)
                theVisualRFmaps(iRGCindex,:,:) = theVisualRFmaps(iRGCindex,:,:) + ...
                    theHartleyPattern * diffResponse(iRGCindex, iFrame);
            end
            
        end  % iFrame


        % Normalize the RFs
        theVisualRFmaps = 1/(2*nStim)*theVisualRFmaps;

        % Compute retinal RFmaps
        marginDegs = min([0.5 0.4*min(theMidgetRGCmosaic.sizeDegs)]);
        spatialSupportSamplesNum = 256;
        retinalRFcenterMaps = theMidgetRGCmosaic.computeRetinalRFcenterMaps(...
            marginDegs, spatialSupportSamplesNum, ...
            'forRGCindices', rgcIndicesOfAnalyzedRFs);

        xLims = theMidgetRGCmosaic.eccentricityDegs(1) + theMidgetRGCmosaic.sizeDegs(1)*0.5*[-1 1];
        yLims = theMidgetRGCmosaic.eccentricityDegs(2) + theMidgetRGCmosaic.sizeDegs(2)*0.5*[-1 1];

        for iRGCindex = 1:numel(retinalRFcenterMaps)
            r = retinalRFcenterMaps{iRGCindex};
            retinalRFcenterMap = r.centerRF;
            retinalSpatialSupportDegsX = r.spatialSupportDegsX;
            retinalSpatialSupportDegsY = r.spatialSupportDegsY;

            cMap = gray(1024);
            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 1700 850], 'Color', [1 1 1]);
            subplot(1,2,1);
            imagesc(retinalSpatialSupportDegsX, retinalSpatialSupportDegsY, retinalRFcenterMap);
            axis 'image'; axis 'xy';
            set(gca, 'XLim', xLims, 'YLim', yLims, 'FontSize', 20);
            title(sprintf('Retinal RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),1), theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),2)));
            
            subplot(1,2,2);
            theRFmap = squeeze(theVisualRFmaps(iRGCindex,:,:));
            theRFmap = theRFmap / max(abs(theRFmap(:)));

            imagesc(visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(1), ...
                    visualRFspatialSupportDegs+theMidgetRGCmosaic.inputConeMosaic.eccentricityDegs(2), ...
                    theRFmap);
            axis 'image'; axis 'xy';
            set(gca, 'XLim', xLims, 'YLim', yLims, 'FontSize', 20, 'Color', cMap(512,:), 'CLim', [-1 1]);
            title(sprintf('Visual RF map \nX,Y = (%2.1f,%2.1f) degs', ...
                theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),1), theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndicesOfAnalyzedRFs(iRGCindex),2)));
            colormap(cMap);
        end % iRGCindex

        % Save data
        fNameRFmaps = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f_RFmaps.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        fNameRFmaps = fullfile(mappedRFsDir,fNameRFmaps);

        save(fNameRFmaps, 'theMidgetRGCmosaic', 'opticsParams', 'visualRFspatialSupportDegs', ...
                          'retinalRFcenterMaps', 'theVisualRFmaps', 'rgcIndicesOfAnalyzedRFs', '-v7.3');
        fprintf('RFmaps for subject %d at position %d of %d saved to %s\n', ...
            subjectRankOrder, iPos, numel(eccXGrid), fNameRFmaps);

    end % iPos

end

function  visualizeComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY, maxVisualizedRFs)

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    rowsNum = 1; %numel(eccX);
    colsNum = 1; %numel(eccY);
    sv = NicePlot.getSubPlotPosVectors(...
       'colsNum', colsNum, ...
       'rowsNum', rowsNum, ...
       'heightMargin',  0.04, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.10, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.02); 

    for iPos = 1:numel(eccXGrid)

        hFig = figure(1);
        set(hFig, 'Position', [10 10 900 900], 'Color', [1 1 1]);

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Load the computed components data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        load(fName, 'thePSFData', 'theMidgetRGCmosaic');
        
        r = floor((iPos-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(iPos-1,colsNum)+1;
        ax = subplot('Position', sv(r,c).v);

        psfData.supportXdegs = thePSFData.psfSupportXdegs;
        psfData.supportYdegs = thePSFData.psfSupportYdegs;
        psfData.data = thePSFData.vLambdaWeightedData;

        theMidgetRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'maxVisualizedRFs', maxVisualizedRFs, ...
            'withSuperimposedPSF', psfData, ...
            'xRange', 0.75, ...
            'yRange', 0.75, ...
            'fontSize', 20);

        fNamePDF = strrep(fName, '.mat', '.pdf');
        NicePlot.exportFigToPDF(fNamePDF, hFig, 300);
    end

end

function modifyOptics(ZernikeDataBase, subjectRankOrder, newZernikeDataBase, newSubjectRankOrder, ...
    newPupilDiamMM, newWavefrontSpatialSamplesNum, eccX, eccY)
    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', [], ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', newPupilDiamMM, ...
        'wavefrontSpatialSamples', newWavefrontSpatialSamplesNum, ...
        'psfUpsampleFactor', 1 ...
        );

    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)
        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];

        % Change position
        opticsParams.positionDegs = mosaicEccDegs;

        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
        load(fName, 'theMidgetRGCmosaic');


        % Generate PSF for new subject
        opticsParams.examinedSubjectRankOrder = newSubjectRankOrder;
        opticsParams.ZernikeDataBase = newZernikeDataBase;

        [thePSFData,~,~,opticsParams] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theMidgetRGCmosaic.inputConeMosaic, []);

        % Save data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            newZernikeDataBase, newSubjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams', '-v7.3');
        fprintf('Data for subject %d at position %d of %d saved to %s\n', newSubjectRankOrder, iPos, numel(eccXGrid), fName);
    end

end


function generateComponents(ZernikeDataBase, subjectRankOrder, eccX, eccY)
    % Struct with the various optics params
    opticsParams = struct(...
        'positionDegs', [], ...           % (x,y) eccentricity for the PSF, in degrees
        'ZernikeDataBase', ZernikeDataBase, ...
        'examinedSubjectRankOrder', subjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
        'analyzedEye', 'right eye', ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'wavefrontSpatialSamples', 701, ...
        'psfUpsampleFactor', 1 ...
        );


    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccXGrid = eccXGrid(:);
    eccYGrid = eccYGrid(:);

    for iPos = 1:numel(eccXGrid)

        % The position analyzed
        mosaicEccDegs = [eccXGrid(iPos) eccYGrid(iPos)];
        radialEccDegs = sqrt(sum(mosaicEccDegs.^2,2));

        if (radialEccDegs == 0)
            mosaicSizeDegs = 0.3;
        elseif (radialEccDegs <= 2)
            mosaicSizeDegs = 0.4;
        elseif (radialEccDegs <= 3)
            mosaicSizeDegs = 0.6;
        elseif (radialEccDegs <= 4)
            mosaicSizeDegs = 0.8;
        elseif (radialEccDegs <= 6)
            mosaicSizeDegs = 1.0;
        elseif (radialEccDegs <= 8)
            mosaicSizeDegs = 1.4;
        elseif (radialEccDegs <= 10)
            mosaicSizeDegs = 1.5;
        elseif (radialEccDegs <= 12)
            mosaicSizeDegs = 1.6;
        elseif (radialEccDegs <= 14)
            mosaicSizeDegs = 1.8;
        elseif (radialEccDegs <= 16)
            mosaicSizeDegs = 2.0;
        elseif (radialEccDegs <= 20)
            mosaicSizeDegs = 2.2;
        else
            mosaicSizeDegs = 2.5;
        end

        % Change position
        opticsParams.positionDegs = mosaicEccDegs;

        % Generate mRGC mosaic
        theMidgetRGCmosaic = midgetRGCMosaic(...
                        'sourceLatticeSizeDegs', 60, ...
                        'whichEye', opticsParams.analyzedEye, ...
                        'eccentricityDegs', mosaicEccDegs, ...
                        'sizeDegs', mosaicSizeDegs*[1 1] ...
                        );

        % Generate PSF
        [thePSFData, ~, ~, opticsParams] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(...
            opticsParams, theMidgetRGCmosaic.inputConeMosaic, []);

        % Save data
        fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
            ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));

        save(fName, 'thePSFData', 'theMidgetRGCmosaic', 'opticsParams', '-v7.3');
        fprintf('Data for position %d of %d saved to %s\n', iPos, numel(eccXGrid), fName);
    end

end
