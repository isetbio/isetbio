function examineMidgetRFcenterSizeVsPSFsize()

    % Eccentricity in temporal retina
    radialEccExamined = [0 1 2 3 4 6 8 12 16 20 24 30];
    for iEcc = 1:numel(radialEccExamined)
        doIt(radialEccExamined(iEcc));
    end

end

function doIt(radialEcc)

    retinaQuadrant = 'temporal meridian';
    whichEye = 'right eye';
    [horizontalEcc, verticalEcc] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaQuadrant, whichEye);
    eccDegs = [horizontalEcc verticalEcc];

    % Optics subject
    opticsDataBase = 'Artal2012';
    opticsSubjectRankOrder = 10;

    theOpticsParams = struct(...
        'positionDegs', eccDegs, ... 
        'ZernikeDataBase', 'Artal2012', ...
        'examinedSubjectRankOrder', 10, ...
        'refractiveErrorDiopters', 0.0, ... 
        'analyzedEye', whichEye, ...
        'subjectRankingEye', 'right eye', ...
        'pupilDiameterMM', 3.0, ...
        'psfUpsampleFactor', 1, ...
        'wavefrontSpatialSamples', 501 ...
        );

    % Compute midget RF center and PSF for this subject at this eccentricity
    destinationRFoverlapRatio = 0.0;
    analyzedCellsNum = 20;
    fName = sprintf('RFcenterPSFanalysis_%2.2fdegs_%2.2fdegs_overlapRatio_%2.2f', horizontalEcc, verticalEcc, destinationRFoverlapRatio);
    [theMidgetRFcenterRcDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeMidgetRFcenterAndPSF(theOpticsParams, destinationRFoverlapRatio, analyzedCellsNum, fName);

    % Save computed radii
    save(sprintf('%s.mat', fName), 'theMidgetRFcenterRcDegs', 'thePSFcharacteristicRadiiDegs', 'exclusiveCenterConesNum', 'opticsDataBase', 'opticsSubjectRankOrder');
end



function [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeMidgetRFcenterAndPSF(theOpticsParams, destinationRFoverlapRatio, analyzedCellsNum, fName)

    % Generate cone mosaic and source lattice struct at the target eccentricity
    targetEccDegs = theOpticsParams.positionDegs;
    targetSizeDegs = 0.5+(sqrt(sum(targetEccDegs.^2,2))+1)*0.4*[0.5 0.25];
    maxVisualizedFOVdegs = round(100 * max(targetSizeDegs) / 6)/100

    [theConeMosaic, sourceLatticeStruct, coneIndicesToBeConnected, thePSFData, theOI] = ...
        generateConeMosaic(targetEccDegs, targetSizeDegs, theOpticsParams);

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
    [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeRFcenters(...
        theMidgetRGCconnectorOBJ, theConeMosaic, thePSFData, ...
        flatTopGaussian, simulateCronerKaplanEstimation, targetEccDegs, maxVisualizedFOVdegs, analyzedCellsNum, fName);

end

function  [retinalRFcenterCharacteristicRadiiDegs, thePSFcharacteristicRadiiDegs, exclusiveCenterConesNum] = analyzeRFcenters(...
    theMidgetRGCconnectorOBJ, theConeMosaic, thePSFData, flatTopGaussian, ...
    simulateCronerKaplanEstimation, targetEccDegs, maxVisualizedFOVdegs, analyzedCellsNum, fName)

    % Compute the center RF sizes for all generated mRGCs
    figNo = 1000;
    hFig = theMidgetRGCconnectorOBJ.visualizeCurrentConnectivity(figNo);
    
    spatialSupportDegs(:,1) = thePSFData.supportXdegs;
    spatialSupportDegs(:,2) = thePSFData.supportYdegs;

    multiStartsNum = 16;
    thePSFData.data = thePSFData.data / max(thePSFData.data(:));
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

    [~,midRow] = min(abs(thePSFData.supportYdegs));
    [~,midCol] = min(abs(thePSFData.supportXdegs));

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
    ax = subplot('Position', subplotPosVectors(1,1).v);
    imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);
    hold(ax, 'on');
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius.x, fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius1.x, fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius2.x, fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
    shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*actualPSFProfileX , yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*actualPSFProfileY, thePSFData.supportYdegs, xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    hold(ax, 'off');
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'CLim', maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {},  'FontSize', 14);
    title(ax, sprintf('v-Lambda weighted PSF\neccentricity (degs): x=%2.1f, y=%2.1f', targetEccDegs(1), targetEccDegs(2)));
    ylabel('space (degs)');
    grid(ax, 'on');


    ax = subplot('Position', subplotPosVectors(1,2).v);
    imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, psfGaussianFit);
    hold(ax, 'on');
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius.x, fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius1.x, fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius2.x, fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
    shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.4*(yLims(2)-yLims(1))*fittedPSFEllipsoidProfileX , yLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    shadedAreaPlot(ax,xLims(1)+0.4*(xLims(2)-xLims(1))*fittedPSFEllipsoidProfileY, thePSFData.supportYdegs, xLims(1), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    hold(ax, 'off');
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'CLim', maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 14);
    title(ax, sprintf('Gaussian fit\nRc (arcmin): %2.2f (%2.2f, %2.2f)', 60*psfCharacteristicRadiusDegs, 60*thePSFcharacteristicRadiiDegs(1), 60*thePSFcharacteristicRadiiDegs(2)));
    grid(ax, 'on');


    ax = subplot('Position', subplotPosVectors(1,3).v);
    imagesc(ax,thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data-psfGaussianFit);
    hold(ax, 'on');
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius.x, fittedPSFEllipsoidOutlineCharacteristicRadius.y, 'b-', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius1.x, fittedPSFEllipsoidOutlineCharacteristicRadius1.y, 'k--', 'LineWidth', 1.0);
    plot(ax, fittedPSFEllipsoidOutlineCharacteristicRadius2.x, fittedPSFEllipsoidOutlineCharacteristicRadius2.y, 'k--', 'LineWidth', 1.0);
    shadedAreaPlot(ax,thePSFData.supportXdegs,yLims(1)+0.1*(yLims(2)-yLims(1)) + 0.4*(yLims(2)-yLims(1))*residualPSFProfileX , yLims(1)+0.1*(yLims(2)-yLims(1)), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    shadedAreaPlot(ax,xLims(1)+0.1*(xLims(2)-xLims(1))+ 0.4*(xLims(2)-xLims(1))*residualPSFProfileY, thePSFData.supportYdegs, xLims(1)+0.1*(xLims(2)-xLims(1)), [0.9 0.9 0.5], [0.5 0.5 0.1], 0.4, 1.0, '--');
    hold(ax, 'off');
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'CLim', maxPSF*[-1 1], 'XLim', xLims, 'YLim', yLims, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 14);
    title(ax, sprintf('residual \n'));
    grid(ax, 'on');
    colormap(brewermap(1024, '*RdBu'));
    drawnow;




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


function [theInputConeMosaic, sourceLatticeStruct, coneIndicesToBeConnected, thePSFData, theOI] = ...
              generateConeMosaic(eccDegs, sizeDegs, theOpticsParams)

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

    % Generate sourceLatticeStruct for the @coneToMidgetRGCConnector
    metaDataStruct.coneTypes = theInputConeMosaic.coneTypes;
    metaDataStruct.coneTypeIDs = [theInputConeMosaic.LCONE_ID theInputConeMosaic.MCONE_ID theInputConeMosaic.SCONE_ID];
    metaDataStruct.coneColors = [theInputConeMosaic.lConeColor; theInputConeMosaic.mConeColor; theInputConeMosaic.sConeColor];

    % Estimate surround radius degs 
    [surroundRadiusDegs, thePSFData, theOI] = estimateSurroundRadiusAtMaxEccentricity(theInputConeMosaic, theOpticsParams);
    metaDataStruct.midgetRGCSurroundRadiusMicronsAtMaxEccentricityGivenOptics = ...
        1e3 * customDegsToMMsConversionFunction(surroundRadiusDegs);

    % Source lattice (i.e. cone mosaic lattice) struct
    sourceLatticeStruct = struct(...
        'name', 'cone RFs', ...
        'DegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
        'MMsToDegsConversionFunction', customMMsToDegsConversionFunction, ...
        'RFpositionsMicrons', theInputConeMosaic.coneRFpositionsMicrons, ...
        'metaData', metaDataStruct ...
        ); 

    coneIndicesToBeConnected = [...
        theInputConeMosaic.lConeIndices(:); ...
        theInputConeMosaic.mConeIndices(:)];

end



function [surroundRadiusDegs, thePSFData, theOI] = estimateSurroundRadiusAtMaxEccentricity(theConeMosaic, theOpticsParams)

    psfWavelengthSupport = [];
    [thePSFData, ~, ~, theOI] = RetinaToVisualFieldTransformer.computeVlambdaWeightedPSF(theOpticsParams, theConeMosaic, psfWavelengthSupport);
    

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
