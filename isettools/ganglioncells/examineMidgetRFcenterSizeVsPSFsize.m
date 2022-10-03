function examineMidgetRFcenterSizeVsPSFsize()

    % Choose retinal quadrant
    retinaQuadrant = 'nasal meridian';
    if (strcmp(retinaQuadrant, 'nasal meridian'))
        radialEccExamined = [1 2 3 4 6 8 12 19 24 30];
    else
        radialEccExamined = [0 1 2 3 4 6 8 12 16 20 24 30];
    end


    % Optics subject
    opticsDataBase = 'Artal2012';

    % When the subject rank order is not empty, we analyze the PSF for that
    % subject and the midget RGC RF center
    opticsSubjectRankOrder = 10;

    % When the subject rank order is empty we analyze the PSFs of all
    % subjects
    %opticsSubjectRankOrder = [];


    % RGC overlap ratio (0 = no overlap)
    destinationRFoverlapRatio = 0.0;


    % If summarizeData is false, we analyze PSF and midget RGC data for the chosen eccentricities
    % If summarizeData is true, we import analyzed PSF and midget RGC data (across the chosen eccentricities) and plot them
    summarizeData = true;
    
    % Do it
%     for iEcc = 1:numel(radialEccExamined)
%         doIt(radialEccExamined(iEcc), summarizeData, opticsDataBase, opticsSubjectRankOrder, retinaQuadrant, destinationRFoverlapRatio);
%     end

    % Plot summarized data
    if (summarizeData)
        dListSingleSubjectPSF = cell(1, numel(radialEccExamined));
        dListAllSubjectPSFs = cell(1, numel(radialEccExamined));
        for iEcc = 1:numel(radialEccExamined)
            dListSingleSubjectPSF{iEcc} = doIt(radialEccExamined(iEcc), summarizeData, opticsDataBase, 10, retinaQuadrant,  destinationRFoverlapRatio);
            dListAllSubjectPSFs{iEcc} = doIt(radialEccExamined(iEcc), summarizeData, opticsDataBase, [], retinaQuadrant, destinationRFoverlapRatio);
        end
    
        plotSummaryData(radialEccExamined, dListSingleSubjectPSF, dListAllSubjectPSFs);
    end
end



function plotSummaryData(radialEccExamined, dListSingleSubjectPSF, dListAllSubjectPSFs)

    exclusiveCenterConesNum = zeros(numel(radialEccExamined),1);
    midgetRGCRFcenterRadiiDegs = zeros(numel(radialEccExamined), 20, 2);
    opticalPSFRadiiDegs = zeros(numel(radialEccExamined),2);
    allSubjectsOpticalPSFRadiiDegs = zeros(numel(radialEccExamined),55,2);

    for iEcc = 1:numel(radialEccExamined)

        % Single subject data - use this for the midget info
        d = dListSingleSubjectPSF{iEcc}
        d.exclusiveCenterConesNum
        % Find cells that were analyzed (20)
        theRcDegs2 = prod(d.midgetRGCRFcenterRadiiDegs,2);
        idx = find(theRcDegs2>0);
        midgetRGCRFcenterRadiiDegs(iEcc,:,:) = d.midgetRGCRFcenterRadiiDegs(idx,:);

        exclusiveCenterConesNum(iEcc,:) = mean(d.exclusiveCenterConesNum);

        % PSF data from the single subject
        opticalPSFRadiiDegs(iEcc,:) = d.opticalPSFRadiiDegs;

        % PSF data from all subjects
        d = dListAllSubjectPSFs{iEcc};
        allSubjectsOpticalPSFRadiiDegs(iEcc,:,:) = d.allSubjectPSFcharacteristicRadiiDegs;
    end

    eccRange = [0-0.5 max(radialEccExamined)+0.5];
    characteristicRadiusRange = [0 7];
    markerSize = 10;

    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 1000 960 960]);
    ax = subplot(2,2,2);
    x = linspace(-1,1,1000);
    Rc = 0.15;
    yGauss = exp(-(x/Rc).^2);
    yGaussPower = exp(-(x/Rc).^4);
    p1 = plot(ax, x, yGauss, 'k-', 'LineWidth', 1.0);
    hold(ax, 'on');
    p2 = plot(ax, x, yGaussPower, 'k-', 'LineWidth', 2, 'Color', [1.0 0.4 0.4]);
    plot(ax, Rc*[-1 1], exp(-1)*[1 1], 'ro-', 'LineWidth', 1.0, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.3 0.3]);
    legend(ax, [p1 p2], {'Gaussian', 'super Gaussian (P=2)'}, ...
        'Location', 'North', 'NumColumns', 2, 'LineWidth', 0.5, 'Color', [0.95 0.95 0.95]);
    set(ax, 'XLim', [-0.5 0.5], 'YLim', [0 1.0], 'XTick', Rc*(-3:1:3), 'XTickLabel', {'-3Rc', '-2Rc', '-Rc', 0, 'Rc', '2Rc', '3Rc'}, ...
        'YTick', 0:0.1:1, 'FontSize', 14);
    axis(ax, 'square'); grid(ax, 'on'); box(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    xlabel(ax,'space')
    ylabel(ax,'sensitivity');

    ax = subplot(2,2,3);
    plot(ax,radialEccExamined, mean(exclusiveCenterConesNum,2), 'ro-', 'LineWidth', 1.5, ...
        'MarkerSize', markerSize, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
    set(ax, 'XLim', eccRange, 'YLim', [0 8], 'XTick', 0:2:30, 'YTick', 0:1:8, 'FontSize', 14);
    axis(ax, 'square'); grid(ax, 'on'); box(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    xlabel(ax,'eccentricity (degs)')
    ylabel(ax,'# of cones in RF center');


    ax = subplot(2,2,4);
    
    allSubjectsMinorCharacteristicRadiiDegs = squeeze(allSubjectsOpticalPSFRadiiDegs(:,:,1));
    allSubjectsMajorCharacteristicRadiiDegs = squeeze(allSubjectsOpticalPSFRadiiDegs(:,:,2));

    shadedAreaBetweenTwoLines(ax,radialEccExamined, ...
        mean(allSubjectsMinorCharacteristicRadiiDegs,2)'*60, ...
        mean(allSubjectsMajorCharacteristicRadiiDegs,2)'*60, ...
        [0.2 1 0.8], [0.1 0.6 0.4], 0.5, 1.0, '-');
    hold(ax, 'on');


%     % Single subject
%     plot(ax,radialEccExamined, max(opticalPSFRadiiDegs,[],2)*60, 'ko-', 'LineWidth', 1.5, ...
%         'MarkerSize', markerSize, 'MarkerFaceColor', [1 1 0.5], 'MarkerEdgeColor', [0.5 0.5 0], 'Color', [0.5 0.5 0]);
%     plot(ax,radialEccExamined, min(opticalPSFRadiiDegs,[],2)*60, 'ko-', 'LineWidth', 1.5, ...
%         'MarkerSize', markerSize, 'MarkerFaceColor', [1 0.7 0.5], 'MarkerEdgeColor', [0.5 0.25 0], 'Color', [0.5 0.25 0]);


    averageRGCsMinorCharacteristicRadiusDegs = mean(midgetRGCRFcenterRadiiDegs(:,:,1),2);
    averageRGCsMajorCharacteristicRadiusDegs = mean(midgetRGCRFcenterRadiiDegs(:,:,2),2);
    averageRGCCharacteristicRadiusDegs = 0.5*(averageRGCsMinorCharacteristicRadiusDegs+averageRGCsMajorCharacteristicRadiusDegs);

    plot(ax, radialEccExamined, averageRGCCharacteristicRadiusDegs*60, 'ro-', 'LineWidth', 1.5, ...
        'MarkerSize', markerSize, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
   

    hold(ax, 'off');
    h = legend(ax,{sprintf('PSF (mean of %d subjects)', size(allSubjectsOpticalPSFRadiiDegs,2)), 'mRGC RF center'}, ...
        'Location', 'NorthWest', 'NumColumns', 1, 'LineWidth', 0.5, 'Color', [0.95 0.95 0.95]);
    
    set(ax, 'XLim', eccRange, 'YLim', characteristicRadiusRange, ...
        'XTick', 0:2:30, 'YTick', 0:1:7, 'FontSize', 14);
    axis(ax, 'square'); grid(ax, 'on'); box(ax, 'off');
    set(ax, 'LineWidth', 1.0, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
    xlabel(ax,'eccentricity (degs)')
    ylabel(ax,'characteristic radius (arcmin)');
end


function summarizedData = doIt(radialEcc, summarizeData, opticsDataBase, opticsSubjectRankOrder, retinaQuadrant, destinationRFoverlapRatio)

    whichEye = 'right eye';
    [horizontalEcc, verticalEcc] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaQuadrant, whichEye);
    eccDegs = [horizontalEcc verticalEcc]


    theOpticsParams = struct(...
        'positionDegs', eccDegs, ... 
        'ZernikeDataBase', opticsDataBase, ...
        'examinedSubjectRankOrder',opticsSubjectRankOrder, ...
        'refractiveErrorDiopters', 0.0, ... 
        'analyzedEye', whichEye, ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'psfUpsampleFactor', 1, ...
        'wavefrontSpatialSamples', 501 ...
        );

    if (isempty(theOpticsParams.examinedSubjectRankOrder))
        fName = sprintf('allSubjectsPSFanalysis_%2.2fdegs_%2.2fdegs_%s', horizontalEcc, verticalEcc, retinaQuadrant);
        if (summarizeData)
            % Load computed radii 
            load(sprintf('%s.mat', fName), 'allSubjectPSFcharacteristicRadiiDegs', 'opticsDataBase', 'opticsSubjectRankOrders');
        else
            [allSubjectPSFcharacteristicRadiiDegs, opticsSubjectRankOrders] = ...
                analyzeAllSubjectPSFs(theOpticsParams, fName);
            % Save computed radii
            save(sprintf('%s.mat', fName), 'allSubjectPSFcharacteristicRadiiDegs', 'opticsDataBase', 'opticsSubjectRankOrders');
        end

        summarizedData = struct();
        summarizedData.allSubjectPSFcharacteristicRadiiDegs = allSubjectPSFcharacteristicRadiiDegs;
        summarizedData.opticsDataBase = opticsDataBase;
        summarizedData.subjectRankOrders = opticsSubjectRankOrders;
 
        return;
    end

    % Compute midget RF center and PSF for this subject at this eccentricity
    analyzedCellsNum = 20;
    fName = sprintf('RFcenterPSFanalysis_%2.2fdegs_%2.2fdegs_overlapRatio_%2.2f_%s', horizontalEcc, verticalEcc, destinationRFoverlapRatio, retinaQuadrant);
    
    if (summarizeData)
        % Load computed radii
        load(sprintf('%s.mat', fName), 'theMidgetRFcenterRcDegs', 'thePSFcharacteristicRadiiDegs', 'exclusiveCenterConesNum', 'opticsDataBase', 'opticsSubjectRankOrder');
    else
        [theMidgetRFcenterRcDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = ...
            analyzeMidgetRFcenterAndPSF(theOpticsParams, destinationRFoverlapRatio, analyzedCellsNum, fName);
        % Save computed radii
        save(sprintf('%s.mat', fName), 'theMidgetRFcenterRcDegs', 'thePSFcharacteristicRadiiDegs', 'exclusiveCenterConesNum', 'opticsDataBase', 'opticsSubjectRankOrder');
    end

    summarizedData = struct();
    summarizedData.midgetRGCRFcenterRadiiDegs = theMidgetRFcenterRcDegs;
    summarizedData.exclusiveCenterConesNum = exclusiveCenterConesNum;
    summarizedData.opticalPSFRadiiDegs = thePSFcharacteristicRadiiDegs;
    summarizedData.opticsDataBase = opticsDataBase;
    summarizedData.opticsSubjectRankOrder = opticsSubjectRankOrder;

end

function [allSubjectPSFcharacteristicRadiiDegs, opticsSubjectRankOrders] = ...
                analyzeAllSubjectPSFs(theOpticsParams, fName)

    opticsSubjectRankOrders = 1:55;

    % Generate cone mosaic and source lattice struct at the target eccentricity
    targetEccDegs = theOpticsParams.positionDegs;
    targetSizeDegs = 0.5+(sqrt(sum(targetEccDegs.^2,2))+1)*0.4*[0.5 0.25];
    maxVisualizedFOVdegs = round(100 * max(targetSizeDegs) / 4)/100;

    

    xLims = maxVisualizedFOVdegs * 0.5*[-1 1];
    yLims = maxVisualizedFOVdegs * 0.5*[-1 1];
    xTicks = -0.5:0.05:0.5;
    yTicks = xTicks;


    videoOBJ = VideoWriter(fName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();


    hFig = figure(54); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1250 450]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 3, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.07, ...
           'leftMargin',     0.07, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);

    allSubjectPSFcharacteristicRadiiDegs = zeros(numel(opticsSubjectRankOrders),2);
    for iSubj = numel(opticsSubjectRankOrders):-1:1

        % Generate PSF for this subject
        theOpticsParams.examinedSubjectRankOrder = opticsSubjectRankOrders(iSubj);

        thePSFData = generateVLambdaWeightedPSF(targetEccDegs, targetSizeDegs, theOpticsParams);

        dFitStruct = fitPSF(thePSFData);

        allSubjectPSFcharacteristicRadiiDegs(iSubj,:) = dFitStruct.thePSFcharacteristicRadiiDegs;

        ax1 = subplot('Position', subplotPosVectors(1,1).v);
        ax2 = subplot('Position', subplotPosVectors(1,2).v);
        ax3 = subplot('Position', subplotPosVectors(1,3).v);
    

        visualizePSFanalysis(ax1, ax2, ax3, thePSFData, dFitStruct, ...
                              xLims, yLims, xTicks, yTicks, theOpticsParams.examinedSubjectRankOrder, targetEccDegs);

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end

    videoOBJ.close();

end


function visualizePSFanalysis(ax1, ax2, ax3, thePSFData, dFitStruct, ...
                              xLims, yLims, xTicks, yTicks, examinedSubjectRankOrder, targetEccDegs)

        ax = ax1;
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);
        hold(ax, 'on');
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*dFitStruct.actualPSFProfileX , yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*dFitStruct.actualPSFProfileY, thePSFData.supportYdegs, xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', dFitStruct.maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {},  'FontSize', 14);
        title(ax, sprintf('v-Lambda weighted PSF (subj rank: %d)\neccentricity (degs): x=%2.1f, y=%2.1f', examinedSubjectRankOrder, targetEccDegs(1), targetEccDegs(2)));
        ylabel('space (degs)');
        grid(ax, 'on');

        ax = ax2;
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, dFitStruct.psfGaussianFit);
        hold(ax, 'on');
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*dFitStruct.fittedPSFEllipsoidProfileX , yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*dFitStruct.fittedPSFEllipsoidProfileY, thePSFData.supportYdegs, xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', dFitStruct.maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 14);
        title(ax, sprintf('Gaussian fit\nRc (arcmin): %2.2f (%2.2f, %2.2f)', ...
            60*dFitStruct.psfCharacteristicRadiusDegs, 60*dFitStruct.thePSFcharacteristicRadiiDegs(1), 60*dFitStruct.thePSFcharacteristicRadiiDegs(2)));
        grid(ax, 'on');

        ax = ax3;
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data-dFitStruct.psfGaussianFit);
        hold(ax, 'on');
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.x, dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.1*(yLims(2)-yLims(1)) + 0.4*(yLims(2)-yLims(1))*dFitStruct.residualPSFProfileX , yLims(1)+0.1*(yLims(2)-yLims(1)), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.1*(xLims(2)-xLims(1))+ 0.4*(xLims(2)-xLims(1))*dFitStruct.residualPSFProfileY, thePSFData.supportYdegs, xLims(1)+0.1*(xLims(2)-xLims(1)), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', dFitStruct.maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 14);
        title(ax, sprintf('residual \n'));
        grid(ax, 'on');
        colormap(brewermap(1024, '*RdBu'));
        drawnow;
end




function [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = ...
    analyzeMidgetRFcenterAndPSF(theOpticsParams, destinationRFoverlapRatio, analyzedCellsNum, fName)

    % Generate cone mosaic and source lattice struct at the target eccentricity
    targetEccDegs = theOpticsParams.positionDegs;
    targetSizeDegs = 0.5+(sqrt(sum(targetEccDegs.^2,2))+1)*0.4*[0.5 0.25];
    maxVisualizedFOVdegs = round(100 * max(targetSizeDegs) / 4)/100;

    [theConeMosaic, sourceLatticeStruct, coneIndicesToBeConnected, thePSFData] = ...
        generateConeMosaicPSFandSourceDestinationLattices(targetEccDegs, targetSizeDegs, theOpticsParams);

    if (isempty(theConeMosaic))
        dFitStruct = fitPSF(thePSFData);

        fprintf('No cones at this eccentricity. Only analyzing PSF data\n');
        retinalRFcenterCharacteristicRadiiDegs = [];
        exclusiveCenterConesNum = 0;
        return;

    end

    % Visualize generated cone mosaic
    theConeMosaic.visualize();

    % Generate the destination lattice struct
    destinationLatticeStruct = generateMRGCMosaicLatticeStruct(...
        sourceLatticeStruct, ...
        theOpticsParams.analyzedEye, ...
        theConeMosaic.customDegsToMMsConversionFunction, ...
        theConeMosaic.customMMsToDegsConversionFunction);

    % [0: minimize chromatic variance 1: minimize spatial variance]
    chromaticSpatialVarianceTradeoff = 1.0;

    % Connect the source and destination lattices
    theMidgetRGCconnectorOBJ = coneToMidgetRGCConnector(...
        sourceLatticeStruct, destinationLatticeStruct, ...
        'optimizationCenter', 'origin', ...
        'chromaticSpatialVarianceTradeoff',chromaticSpatialVarianceTradeoff, ...
        'coneIndicesToBeConnected', coneIndicesToBeConnected, ...
        'visualizeConnectivityAtIntermediateStages', ~true, ...
        'generateProgressVideo', ~true);

    % Diverge input cones to multiple nearby midget RGCs (RF overlap)
    theMidgetRGCconnectorOBJ.divergeSourceRFsToNearbyDestinationRFs('destinationRFoverlapRatio', destinationRFoverlapRatio);

    flatTopGaussian = true;
    simulateCronerKaplanEstimation = ~true;
    [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeRFcentersAndPSF(...
        theMidgetRGCconnectorOBJ, theConeMosaic, thePSFData, theOpticsParams.examinedSubjectRankOrder, ...
        flatTopGaussian, simulateCronerKaplanEstimation, targetEccDegs, maxVisualizedFOVdegs, analyzedCellsNum, fName);

end


function dFitStruct = fitPSF(thePSFData)

    spatialSupportDegs(:,1) = thePSFData.supportXdegs;
    spatialSupportDegs(:,2) = thePSFData.supportYdegs;

    multiStartsNum = 16;
    theFittedGaussianPSF = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportDegs(:,1), spatialSupportDegs(:,2), thePSFData.data, ...
            'flatTopGaussian', false, ...
            'forcedOrientationDegs', [], ...
            'forcedCentroidXYpos', [], ...
            'globalSearch', true, ...
            'multiStartsNum', multiStartsNum);
    
    % The sqrt(product) of the 2 radii
    thePSFcharacteristicRadiiDegs = theFittedGaussianPSF.characteristicRadii;
    psfCharacteristicRadiusDegs = sqrt(prod(thePSFcharacteristicRadiiDegs));
    
    % The fitted Gaussian ellipsoid fit
    psfGaussianFit = theFittedGaussianPSF.ellipsoidMap;

    psfGaussianCenter = theFittedGaussianPSF.xyCenter;

    maxPSF = max([max(thePSFData.data(:)) max(psfGaussianFit(:))]);

    [~,idx] = max(psfGaussianFit(:));
    [midPSFRow, midPSFCol] = ind2sub(size(psfGaussianFit), idx);


    fittedPSFEllipsoidProfileX = squeeze(psfGaussianFit(midPSFRow,:))/max(psfGaussianFit(:));
    fittedPSFEllipsoidProfileY = squeeze(psfGaussianFit(:,midPSFCol))/max(psfGaussianFit(:));
    actualPSFProfileX = squeeze(thePSFData.data(midPSFRow,:))/max(psfGaussianFit(:));
    actualPSFProfileY = squeeze(thePSFData.data(:,midPSFCol))/max(psfGaussianFit(:));
    residualPSFProfileX = squeeze(thePSFData.data(midPSFRow,:)-psfGaussianFit(midPSFRow,:))/max(psfGaussianFit(:));
    residualPSFProfileY = squeeze(thePSFData.data(:,midPSFCol)-psfGaussianFit(:,midPSFCol))/max(psfGaussianFit(:));

    iAngles = 0:10:360;
    fittedPSFEllipsoidOutlineCharacteristicRadius.x = psfGaussianCenter(1) + psfCharacteristicRadiusDegs * cosd(iAngles);
    fittedPSFEllipsoidOutlineCharacteristicRadius.y = psfGaussianCenter(2) + psfCharacteristicRadiusDegs * sind(iAngles);

    fittedPSFEllipsoidOutlineCharacteristicRadius1.x = psfGaussianCenter(1) + thePSFcharacteristicRadiiDegs(1) * cosd(iAngles);
    fittedPSFEllipsoidOutlineCharacteristicRadius1.y = psfGaussianCenter(2) + thePSFcharacteristicRadiiDegs(1) * sind(iAngles);

    fittedPSFEllipsoidOutlineCharacteristicRadius2.x = psfGaussianCenter(1) + thePSFcharacteristicRadiiDegs(2) * cosd(iAngles);
    fittedPSFEllipsoidOutlineCharacteristicRadius2.y = psfGaussianCenter(2) + thePSFcharacteristicRadiiDegs(2) * sind(iAngles);

    % Gather everything in a struct
    dFitStruct.thePSFcharacteristicRadiiDegs = thePSFcharacteristicRadiiDegs;
    dFitStruct.psfCharacteristicRadiusDegs = psfCharacteristicRadiusDegs;
    dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius = fittedPSFEllipsoidOutlineCharacteristicRadius;
    dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius1 = fittedPSFEllipsoidOutlineCharacteristicRadius1;
    dFitStruct.fittedPSFEllipsoidOutlineCharacteristicRadius2 = fittedPSFEllipsoidOutlineCharacteristicRadius2;
    dFitStruct.fittedPSFEllipsoidProfileX = fittedPSFEllipsoidProfileX;
    dFitStruct.fittedPSFEllipsoidProfileY = fittedPSFEllipsoidProfileY;
    dFitStruct.actualPSFProfileX = actualPSFProfileX;
    dFitStruct.actualPSFProfileY = actualPSFProfileY;
    dFitStruct.residualPSFProfileX = residualPSFProfileX;
    dFitStruct.residualPSFProfileY = residualPSFProfileY;
    dFitStruct.psfGaussianFit =  psfGaussianFit;
    dFitStruct.maxPSF = maxPSF;

end


function  [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeRFcentersAndPSF(...
    theMidgetRGCconnectorOBJ, theConeMosaic, thePSFData, examinedSubjectRankOrder, flatTopGaussian, ...
    simulateCronerKaplanEstimation, targetEccDegs, maxVisualizedFOVdegs, analyzedCellsNum, fName)

    % Compute the center RF sizes for all generated mRGCs
    figNo = 1000;
    hFig = theMidgetRGCconnectorOBJ.visualizeCurrentConnectivity(figNo);
    
    % Fit the PSF data
    dFitStruct = fitPSF(thePSFData);
    thePSFcharacteristicRadiiDegs = dFitStruct.thePSFcharacteristicRadiiDegs;

    
    spatialSupportDegs(:,1) = thePSFData.supportXdegs;
    spatialSupportDegs(:,2) = thePSFData.supportYdegs;

    [~,midRow] = min(abs(thePSFData.supportYdegs));
    [~,midCol] = min(abs(thePSFData.supportXdegs));


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.07, ...
           'leftMargin',     0.07, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);

    xLims = maxVisualizedFOVdegs * 0.5*[-1 1];
    yLims = maxVisualizedFOVdegs * 0.5*[-1 1];
    xTicks = -0.5:0.05:0.5;
    yTicks = xTicks;


    videoOBJ = VideoWriter(fName, 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();


    hFig = figure(55); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1250 900]);

    % The PSF
    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    ax2 = subplot('Position', subplotPosVectors(1,2).v);
    ax3 = subplot('Position', subplotPosVectors(1,3).v);
    visualizePSFanalysis(ax1, ax2, ax3, thePSFData, dFitStruct, ...
                              xLims, yLims, xTicks, yTicks, examinedSubjectRankOrder, targetEccDegs);



    theCentroids = theMidgetRGCconnectorOBJ.destinationRFcentroidsFromInputs;
    theMeanCentroid = mean(theCentroids,1);
    [~,idx] = sort(sum((bsxfun(@minus, theCentroids, theMeanCentroid)).^2,2));

    
    visualizedRGCindices = idx(1:min([analyzedCellsNum numel(idx)]));

    retinalRFcenterCharacteristicRadiiDegs = zeros(numel(visualizedRGCindices),2);
    exclusiveCenterConesNum = zeros(numel(visualizedRGCindices),1);

    fprintf('Analyzing the central %d of a total of %d RGCs\n', numel(visualizedRGCindices), size(theMidgetRGCconnectorOBJ.connectivityMatrix,2));
    


    % Analyze the central-most RFs
    for iSortedRGCindex = 1:numel(visualizedRGCindices)

        mRGCRFindex = visualizedRGCindices(iSortedRGCindex);

        theConeIndices = find(abs(squeeze(theMidgetRGCconnectorOBJ.connectivityMatrix(:, mRGCRFindex)))>0);
        theConeWeights = abs(full(theMidgetRGCconnectorOBJ.connectivityMatrix(theConeIndices, mRGCRFindex)));

        theExclusiveConeIndices = find(squeeze(theMidgetRGCconnectorOBJ.connectivityMatrix(:, mRGCRFindex))>0);
        exclusiveCenterConesNum(iSortedRGCindex) = numel(theExclusiveConeIndices);

        [retinalRFcenterConeMap, retinalRFcenterConeMapGaussianFit, retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,:), retinalRFcenterDegs, ...
            flatTopExponents, meanConePositionDegs] = ...
            computeAnatomicalRF(theConeMosaic, theConeWeights, theConeIndices, spatialSupportDegs, simulateCronerKaplanEstimation, flatTopGaussian);
        
        retinalRFcenterCharacteristicRadiusDegs = sqrt(prod(retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,:),2));

        maxRF = max([max(retinalRFcenterConeMap(:)) max(retinalRFcenterConeMapGaussianFit(:))]);

        iAngles = 0:10:360;
        fittedGaussianEllipsoidOutlineCharacteristicRadius.x = ...
            retinalRFcenterDegs(1)-meanConePositionDegs(1) + retinalRFcenterCharacteristicRadiusDegs  * cosd(iAngles);
        fittedGaussianEllipsoidOutlineCharacteristicRadius.y = ...
            retinalRFcenterDegs(2)-meanConePositionDegs(2) + retinalRFcenterCharacteristicRadiusDegs  * sind(iAngles);
        
        fittedGaussianEllipsoidOutlineCharacteristicRadius1.x = ...
            retinalRFcenterDegs(1)-meanConePositionDegs(1) + retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,1) * cosd(iAngles);
        fittedGaussianEllipsoidOutlineCharacteristicRadius1.y = ...
            retinalRFcenterDegs(2)-meanConePositionDegs(2) + retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,1) * sind(iAngles);

        fittedGaussianEllipsoidOutlineCharacteristicRadius2.x = ...
            retinalRFcenterDegs(1)-meanConePositionDegs(1) + retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,2) * cosd(iAngles);
        fittedGaussianEllipsoidOutlineCharacteristicRadius2.y = ...
            retinalRFcenterDegs(2)-meanConePositionDegs(2) + retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,2) * sind(iAngles);


        fittedGaussianEllipsoidProfileX = squeeze(retinalRFcenterConeMapGaussianFit(midRow,:))/max(retinalRFcenterConeMapGaussianFit(:));
        fittedGaussianEllipsoidProfileY = squeeze(retinalRFcenterConeMapGaussianFit(:,midCol))/max(retinalRFcenterConeMapGaussianFit(:));


        ax = subplot('Position', subplotPosVectors(2,1).v);
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, retinalRFcenterConeMap);
        hold(ax, 'on');
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius.x, fittedGaussianEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius1.x, fittedGaussianEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius2.x, fittedGaussianEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*fittedGaussianEllipsoidProfileX, yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*fittedGaussianEllipsoidProfileY, thePSFData.supportYdegs,xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', maxRF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'FontSize', 14);
        title(ax, sprintf('retinal RF center\neccentricity (degs): x=%2.2f, y=%2.2f (%d of %d)', ...
            meanConePositionDegs(1), meanConePositionDegs(2), ...
            iSortedRGCindex, numel(visualizedRGCindices)));
        xlabel('space (degs)');
        ylabel('space (degs)');
        grid(ax, 'on');

        ax = subplot('Position', subplotPosVectors(2,2).v);
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, retinalRFcenterConeMapGaussianFit);
        hold(ax, 'on');
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius.x, fittedGaussianEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius1.x, fittedGaussianEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius2.x, fittedGaussianEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*fittedGaussianEllipsoidProfileX, yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*fittedGaussianEllipsoidProfileY, thePSFData.supportYdegs,xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', max(retinalRFcenterConeMapGaussianFit(:))*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', {}, 'FontSize', 14);
        if (flatTopGaussian)
            title(ax, sprintf('flat-top Gaussian ellipsoid fit (exponents: %2.2f, %2.2f)\nRc (arcmin): %2.2f (%2.2f, %2.2f)', 2*flatTopExponents(1), 2*flatTopExponents(2), ...
                60*retinalRFcenterCharacteristicRadiusDegs, 60*retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,1), 60*retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,2)));
        else
            title(ax, sprintf('Gaussian ellipsoid fit\nRc (arcmin): %2.2f (%2.2f, %2.2f)', ...
                60*retinalRFcenterCharacteristicRadiusDegs, 60*retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,1), 60*retinalRFcenterCharacteristicRadiiDegs(mRGCRFindex,2)));
        end
        xlabel('space (degs)');
        grid(ax, 'on');

        ax = subplot('Position', subplotPosVectors(2,3).v);
        imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, retinalRFcenterConeMap-retinalRFcenterConeMapGaussianFit);
        hold(ax, 'on');
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius.x, fittedGaussianEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius1.x, fittedGaussianEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
        plot(ax, fittedGaussianEllipsoidOutlineCharacteristicRadius2.x, fittedGaussianEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
        shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*fittedGaussianEllipsoidProfileX, yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*fittedGaussianEllipsoidProfileY, thePSFData.supportYdegs,xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'CLim', maxRF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', {}, 'FontSize', 14);
        title(ax, sprintf('residual\n'));
        grid(ax, 'on');
        xlabel('space (degs)');

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end

    videoOBJ.close();
end

function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end


function [retinalRFcenterConeMap, retinalRFcenterConeMapGaussianFit, retinalRFcenterCharacteristicRadiiDegs, retinalRFcenter, flatTopExponents, meanConePositionDegs] = ...
    computeAnatomicalRF(theConeMosaic, theConeWeights, theConeIndices, spatialSupportDegs, simulateCronerKaplanEstimation, flatTopGaussian)

    % Compute the retinal RF center cone map
    theConeCharacteristicRadiiDegs = theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * theConeMosaic.coneApertureDiametersDegs(theConeIndices);
    theConePositionsDegs = theConeMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    spatialSupportDegs = bsxfun(@plus, spatialSupportDegs, meanConePositionDegs);
    retinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConePositionsDegs, theConeWeights, spatialSupportDegs);

    if (simulateCronerKaplanEstimation)
        % Fit a 1D Gaussian line weighting function to the 1D profile 
        % (integration along the Y-dimension of the 2D visually projected
        % cone aperture map)
        retinalRFcenterConeMapProfile = sum(retinalRFcenterConeMap,1);
        retinalRFcenterConeMapProfile = retinalRFcenterConeMapProfile / max(retinalRFcenterConeMapProfile(:)) * max(retinalRFcenterConeMap(:));
        theFittedGaussianLineWeightingFunction = RetinaToVisualFieldTransformer.fitGaussianLineWeightingFunction(...
            squeeze(spatialSupportDegs(:,1)), retinalRFcenterConeMapProfile(:));

        % Return the characteristic radius in degrees
        retinalRFcenterCharacteristicRadiiDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;
        
        % Make a 2D circularly-symmetric version of the fitted profile
        retinalRFcenterConeMapGaussianFitProfile = theFittedGaussianLineWeightingFunction.profile(:);
        retinalRFcenterConeMapGaussianFit = retinalRFcenterConeMapGaussianFitProfile * (retinalRFcenterConeMapGaussianFitProfile)';

        flatTopExponents = [];
    else
        rangeForEllipseRcYRcXratio = [1.0/1.4 1.4];
        forcedOrientationDegs = [];
        globalSearch = true;
        multiStartsNum = 16;
        forcedCentroidXYpos = meanConePositionDegs;

        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportDegs(:,1), spatialSupportDegs(:,2), retinalRFcenterConeMap, ...
            'flatTopGaussian', flatTopGaussian, ...
            'forcedEllipseRcYRcXratio', [], ...
            'rangeForEllipseRcYRcXratio', rangeForEllipseRcYRcXratio, ...
            'forcedOrientationDegs', forcedOrientationDegs, ...
            'forcedCentroidXYpos', forcedCentroidXYpos, ...
            'globalSearch', globalSearch, ...
            'multiStartsNum', multiStartsNum);
    
        % The characteristic radii
        retinalRFcenterCharacteristicRadiiDegs = theFittedGaussian.characteristicRadii;
    
        % The fitted Gaussian ellipsoid fit
        retinalRFcenterConeMapGaussianFit = theFittedGaussian.ellipsoidMap;

        % The center coords of the fitted Gaussian ellipsoid
        retinalRFcenter = theFittedGaussian.xyCenter;

        % The exponents
        flatTopExponents = theFittedGaussian.flatTopExponents;
    end


end



function destinationLatticeStruct = generateMRGCMosaicLatticeStruct(sourceLatticeStruct, whichEye, customDegsToMMsConversionFunction, customMMsToDegsConversionFunction)
    % Destination lattice struct

    sourceLatticeSizeDegs = 60;

    % Import mRGC RF positions (destination) 
    mRGCRFposMicrons = retinalattice.import.finalMRGCPositions(...
                 sourceLatticeSizeDegs, ...
                 mean(sourceLatticeStruct.RFpositionsMicrons,1), ... 
                 max(max(sourceLatticeStruct.RFpositionsMicrons,[], 1)-min(sourceLatticeStruct.RFpositionsMicrons,[], 1)), ...
                 whichEye, ...
                 customDegsToMMsConversionFunction);
    
    destinationLatticeStruct = struct(...
        'name', 'mRGC RFs', ...
        'DegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', mRGCRFposMicrons ...
        );
end


function thePSFData = generateVLambdaWeightedPSF(eccDegs, sizeDegs, theOpticsParams)
    theInputConeMosaic = generateConeMosaic(eccDegs, sizeDegs, theOpticsParams);

    psfWavelengthSupport = [];
    [thePSFData, ~, ~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(theOpticsParams, theInputConeMosaic, psfWavelengthSupport);
    thePSFData.data = thePSFData.data / max(thePSFData.data(:));
end


function theInputConeMosaic = generateConeMosaic(eccDegs, sizeDegs, theOpticsParams)
     % Set cone aperture modifiers
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');

    sourceLatticeSizeDegs = 60;
    customDegsToMMsConversionFunction = @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x);
    customMMsToDegsConversionFunction = @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x);

    % Generate the input cone mosaic
    theInputConeMosaic = cMosaic(...
       'sourceLatticeSizeDegs', sourceLatticeSizeDegs, ...
       'eccentricityDegs', eccDegs, ...
       'sizeDegs', sizeDegs, ...
       'whichEye', theOpticsParams.analyzedEye, ...
       'coneDensities', [0.6 0.3 0.1], ...
       'overlappingConeFractionForElimination', 0.5, ...
       'rodIntrusionAdjustedConeAperture', true, ...
       'coneApertureModifiers', coneApertureModifiers, ...
       'customDegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
       'customMMsToDegsConversionFunction', customMMsToDegsConversionFunction);

end


function [theInputConeMosaic, sourceLatticeStruct, coneIndicesToBeConnected, thePSFData] = ...
              generateConeMosaicPSFandSourceDestinationLattices(eccDegs, sizeDegs, theOpticsParams)

    % Generate the input cone mosaic
    theInputConeMosaic = generateConeMosaic(eccDegs, sizeDegs, theOpticsParams);

    
    psfWavelengthSupport = [];
    [thePSFData, ~, ~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(theOpticsParams, theInputConeMosaic, psfWavelengthSupport);
    thePSFData.data = thePSFData.data / max(thePSFData.data(:));

    if (numel(theInputConeMosaic.coneTypes) == 0)
        theInputConeMosaic = [];
        sourceLatticeStruct = [];
        coneIndicesToBeConnected = [];
        return;
    end

    % Generate sourceLatticeStruct for the @coneToMidgetRGCConnector
    metaDataStruct.coneTypes = theInputConeMosaic.coneTypes;
    metaDataStruct.coneTypeIDs = [theInputConeMosaic.LCONE_ID theInputConeMosaic.MCONE_ID theInputConeMosaic.SCONE_ID];
    metaDataStruct.coneColors = [theInputConeMosaic.lConeColor; theInputConeMosaic.mConeColor; theInputConeMosaic.sConeColor];

    % Estimate surround radius degs 
    surroundRadiusDegs = estimateSurroundRadiusAtMaxEccentricity(theInputConeMosaic, thePSFData);
    metaDataStruct.midgetRGCSurroundRadiusMicronsAtMaxEccentricityGivenOptics = ...
        1e3 * theInputConeMosaic.customDegsToMMsConversionFunction(surroundRadiusDegs);

    % Source lattice (i.e. cone mosaic lattice) struct
    sourceLatticeStruct = struct(...
        'name', 'cone RFs', ...
        'DegsToMMsConversionFunction', theInputConeMosaic.customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', theInputConeMosaic.customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', theInputConeMosaic.coneRFpositionsMicrons, ...
        'metaData', metaDataStruct ...
        ); 

    coneIndicesToBeConnected = [...
        theInputConeMosaic.lConeIndices(:); ...
        theInputConeMosaic.mConeIndices(:)];

end



function surroundRadiusDegs = estimateSurroundRadiusAtMaxEccentricity(theConeMosaic, thePSFData)

    
    % Estimate mean anatomical cone aperture from the 6 closest (to the
    % cone with the max position.
    [~,idx] = MosaicConnector.pdist2(theConeMosaic.coneRFpositionsDegs, [], ...
        'fromPosition', 'maxAbsPosition', ...
        'smallest', 6 ...
        );

    meanConeApertureDegs = mean(theConeMosaic.coneApertureDiametersDegs(idx));
    
    anatomicalConeCharacteristicRadiusDegs = theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor  * meanConeApertureDegs;

    % Do not simulate Croner/Kaplan estimation if we want to fit a Gaussian
    % point spread function (2D). The Croner/Kaplan estimation fits a
    % Gaussian line spread function (1D) 
    simulateCronerKaplanEstimation = false;

    hFig = figure(1); clf;
    visualConeCharacteristicRadiusDegs = RetinaToVisualFieldTransformer.analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, simulateCronerKaplanEstimation, hFig);

    RsRcRatio = 6.7;
    characteristicRadiiLimit = 2;
    surroundRadiusDegs = characteristicRadiiLimit * (visualConeCharacteristicRadiusDegs*RsRcRatio);

end
