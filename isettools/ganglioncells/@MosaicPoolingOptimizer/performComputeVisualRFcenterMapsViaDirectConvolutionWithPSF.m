function performComputeVisualRFcenterMapsViaDirectConvolutionWithPSF()
    
    horizontalEccDegs = [-15 -10 -7 -4 -2 0 2 4 7 10 15];
    verticalEccDegs = [-10 -7 -4 -2 0 2 4 7 10];
    

    horizontalEccDegs = [-8 0 8];
    verticalEccDegs = [-8 0 8];

    [X,Y] = meshgrid(horizontalEccDegs, verticalEccDegs);
    XX = X(:);
    YY = Y(:);

    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';
    
    for iPos = 1:numel(XX)
        [iRow, iCol] = ind2sub(size(X), iPos);

        % Set up figures
        hFig1 = figure(2000); clf;
        ff = MSreadyPlot.figureFormat('1x1 small no labels');
        theAxes = MSreadyPlot.generateAxes(hFig1,ff);
        axRetinal = theAxes{1,1};
        set(hFig1, 'Color', [1 1 1]);

        hFig2 = figure(2001); clf;
        ff = MSreadyPlot.figureFormat('1x1 small no labels');
        theAxes = MSreadyPlot.generateAxes(hFig2,ff);
        axVisual = theAxes{1,1};
        set(hFig2, 'Color', [1 1 1]);

        hFig3 = figure(2002); clf;
        ff = MSreadyPlot.figureFormat('1x1 small no labels');
        theAxes = MSreadyPlot.generateAxes(hFig3,ff);
        axPSF = theAxes{1,1};
        set(hFig3, 'Color', [1 1 1]);

        % Run it !
        runAtEccentricity(axRetinal, axVisual, axPSF, ff, [horizontalEccDegs(iCol) verticalEccDegs(iRow)]);

        % Export PDFs
        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('retinalRFs_%d_%d.pdf', horizontalEccDegs(iCol), verticalEccDegs(iRow)));
        NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig1, 300);

        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('visualRFs_%d_%d.pdf', horizontalEccDegs(iCol), verticalEccDegs(iRow)));
        NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig2, 300);

        pdfFileNameForPLOS = fullfile(rawFiguresRoot, sprintf('PSFs_%d_%d.pdf', horizontalEccDegs(iCol), verticalEccDegs(iRow)));
        NicePlot.exportFigToPDF(pdfFileNameForPLOS, hFig3, 300);
    end

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);
end

function runAtEccentricity(axRetinal, axVisual, axPSF, ff, eccDegs)
    mosaicParams.eccDegs = eccDegs;
    mosaicParams.sizeDegs = [1 1];
    mosaicParams.chromaticSpatialVarianceTradeoff = 1.0;
    mosaicParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = MosaicConnector.maxSourceInputsToConsiderTransferToNearbyDestinationRF;
    theMidgetRGCMosaic = generateTheCenterOnlyConnectedMRGCmosaic(mosaicParams);


    % Optics params
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;

    % Change subject rank order to see other subjects
    %opticsParams.examinedSubjectRankOrder = 

    % 'Native optics', at mosaic's eccentricity
    opticsParams.positionDegs = [];

    % PSF size and resolution
    opticsParams.wavefrontSpatialSamples = 501;
    mosaicRadialEcc = sqrt(sum(mosaicParams.eccDegs(:).^2));
    if (mosaicRadialEcc < 0.3)
        opticsParams.psfUpsampleFactor = 3;
    else
        opticsParams.psfUpsampleFactor = 1;
    end
    

    % Generate the Vlambda weighted psfData
    thePSFData = MosaicPoolingOptimizer.generateVlambdaWeightedPSFData(theMidgetRGCMosaic, opticsParams);

    dd = bsxfun(@minus,theMidgetRGCMosaic.rgcRFpositionsDegs, mosaicParams.eccDegs);
    dd = sqrt(sum(dd.^2,2));
    [~,idx] = sort(dd, 'ascend');

    % Examine the 16 RGCs closest to the mosaic's center
    maxNeuronsVisualized = 19;
    visualizedSizeDegs = 0.4;

    visualizedRFsNum = min([maxNeuronsVisualized numel(idx)]);
    examinedRGCindices = idx(1:visualizedRFsNum);

    generateFigures = false;
    [retinalRFcontourData, visualRFcontourData, xSupportDegs, ySupportDegs] = ...
        analyzeRFmaps(theMidgetRGCMosaic, examinedRGCindices, thePSFData, ...
        generateFigures);

   
    for iRGC = 1:numel(retinalRFcontourData)
        hold(axRetinal, 'on');

        if (~isempty(retinalRFcontourData{iRGC}))
            S = retinalRFcontourData{iRGC};
            S.FaceVertexCData = [0.85 0.8 0.35];
            S.FaceColor = 'flat';
            S.EdgeColor = [0 0 0];
            S.FaceAlpha = 0.3;
            S.LineWidth = 1.0;
            patch(S, 'Parent', axRetinal)
        else
            fprintf(2,'Empty contour. No plotting\n')
        end


        hold(axVisual, 'on');
        if (~isempty(visualRFcontourData{iRGC}))
            S = visualRFcontourData{iRGC};
            S.FaceVertexCData = [0.05 0.35 0.85];
            S.FaceColor = 'flat';
            S.EdgeColor = [0 0 0];
            S.FaceAlpha = 0.3;
            S.LineWidth = 1.0;
            patch(S, 'Parent', axVisual);
        else
            fprintf(2,'Empty contour. No plotting\n')
        end
    end % iRGC

    xlabel(axRetinal, '');
    ylabel(axRetinal, '');

    axis(axRetinal, 'equal'); axis(axRetinal, 'xy');
    set(axRetinal, 'XLim', mosaicParams.eccDegs(1) + visualizedSizeDegs*0.5*[-1 1], ... 
                   'YLim', mosaicParams.eccDegs(2) + visualizedSizeDegs*0.5*[-1 1]);
    set(axRetinal, 'XTick', [], ...
                   'YTick', []);
    box(axRetinal, 'on');
    % axis color and width
    set(axRetinal, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);


    xlabel(axVisual, '');
    ylabel(axVisual, '');
    axis(axVisual, 'equal'); axis(axVisual, 'xy');
    set(axVisual, 'XLim', mosaicParams.eccDegs(1) + visualizedSizeDegs*0.5*[-1 1], ... 
                  'YLim', mosaicParams.eccDegs(2) + visualizedSizeDegs*0.5*[-1 1]);
    set(axVisual, 'XTick', [], ...
                   'YTick', []);
    box(axVisual, 'on');


    cmap = brewermap(1024,'greys');
    alpha = 0.5;
    contourLineColor = [0.0 0.0 0.0];
    cMosaic.semiTransparentContourPlot(axPSF, xSupportDegs,ySupportDegs, ...
        thePSFData.vLambdaWeightedPSF/max(thePSFData.vLambdaWeightedPSF(:)), 0.05:0.2:0.95, cmap, alpha, contourLineColor);
    
    axis(axPSF, 'image');  axis(axPSF, 'xy');
    colormap(axPSF, brewermap(1024, '*greys'));
    
    xlabel(axPSF, '');
    ylabel(axPSF, '');
    set(axPSF, 'XLim', mosaicParams.eccDegs(1) + visualizedSizeDegs*0.5*[-1 1], ... 
                  'YLim', mosaicParams.eccDegs(2) + visualizedSizeDegs*0.5*[-1 1]);
    set(axPSF, 'XTick', [], ...
                   'YTick', []);
    box(axPSF, 'on');
    % axis color and width
    set(axPSF, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
end


function [retinalRFcontourData, visualRFcontourData, xSupportDegs, ySupportDegs] = ...
    analyzeRFmaps(theMidgetRGCMosaic, examinedRGCindices, thePSFData, generateFigures)

    centerSubregionContourSamples = 64;

    xSupportDegs = theMidgetRGCMosaic.eccentricityDegs(1)+thePSFData.supportX/60;
    ySupportDegs = theMidgetRGCMosaic.eccentricityDegs(2)+thePSFData.supportY/60;

    if (generateFigures)
        figure(1); clf;
    end

    retinalRFcontourData = cell(1, numel(examinedRGCindices));
    visualRFcontourData = cell(1, numel(examinedRGCindices));

    % Find theStimulusPixelSizeDegs
    theStimulusPixelSizeDegs = zeros(1, numel(examinedRGCindices));
    for iRGC = 1:numel(examinedRGCindices)
        theRGCindex = examinedRGCindices(iRGC);
        connectivityVector = full(squeeze(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        centerConeIndices = find(connectivityVector> 0.001);
        coneHalfSpacings = theMidgetRGCMosaic.inputConeMosaic.coneRFspacingsDegs(centerConeIndices);
        theStimulusPixelSizeDegs(iRGC) = mean(coneHalfSpacings);
    end
    theStimulusPixelSizeDegs = 0.25 * mean(theStimulusPixelSizeDegs);

    % Compute stimulus pixel image
    thePSFsampleSizeDegs = xSupportDegs(2)-xSupportDegs(1);
    stimulusPixelSamples = max([1 round(theStimulusPixelSizeDegs/thePSFsampleSizeDegs) ]);
    theStimulusPixelImage = ones(stimulusPixelSamples, stimulusPixelSamples);

    for iRGC = 1:numel(examinedRGCindices)
        theRGCindex = examinedRGCindices(iRGC);

        connectivityVector = full(squeeze(theMidgetRGCMosaic.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        centerConeIndices = find(connectivityVector> 0.001);
        poolingWeights = connectivityVector(centerConeIndices);
        conePos = theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,:);
        coneRc = ...
                theMidgetRGCMosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                theMidgetRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(centerConeIndices);

        % Convolve PSF with the image of the RF mapping stimulus image
        theOpticalStimulusPSF = conv2(thePSFData.vLambdaWeightedPSF, theStimulusPixelImage,'same');

        [~,retinalRFmap2D] = mRGCMosaic.subregionOutlineContourFromPooledCones(...
          conePos, coneRc, poolingWeights, ...
          xSupportDegs, ySupportDegs, centerSubregionContourSamples);
    
        theContourData = mRGCMosaic.subregionEllipseFromPooledCones(...
            conePos, coneRc, poolingWeights, ...
            xSupportDegs, ySupportDegs, [], centerSubregionContourSamples);
        
        if (isempty(theContourData))
            retinalRFcontourData{iRGC} = [];
            visualRFcontourData{iRGC} = [];
            continue;
        end


        S = [];
        if (iscell(theContourData))
            S.Vertices = theContourData{1}.vertices;
            S.Faces = theContourData{1}.faces;
        else
            S.Vertices = theContourData.vertices;
            S.Faces = theContourData.faces;
        end
        retinalRFcontourData{iRGC} = S;


        % Compute visual image of retinal RF map
        theVisualRFmap = conv2(retinalRFmap2D, theOpticalStimulusPSF, 'same');


        % Fit an ellipse to theVisualRFmap and draw a contour at its
        % characteristic radius
        zLevel = exp(-1);
        theContourData = mRGCMosaic.ellipseContourFromSubregionRFmap(...
            xSupportDegs, ySupportDegs, theVisualRFmap, zLevel, centerSubregionContourSamples);

        if (isempty(theContourData))
            retinalRFcontourData{iRGC} = [];
            visualRFcontourData{iRGC} = [];
            continue;
        end

        S = [];
        if (iscell(theContourData))
            S.Vertices = theContourData{1}.vertices;
            S.Faces = theContourData{1}.faces;
        else
            S.Vertices = theContourData.vertices;
            S.Faces = theContourData.faces;
        end
        visualRFcontourData{iRGC} = S;


        if (generateFigures)
            ax = subplot(2,3,1);
            imagesc(ax,xSupportDegs,ySupportDegs, retinalRFmap2D);
            colormap(ax, gray(1024))
            hold (ax, 'on');

            S = retinalRFcontourData{iRGC};
            S.FaceVertexCData = [0.5 0.5 0.5];
            S.FaceColor = 'flat';
            S.EdgeColor = [1 0 0];
            S.FaceAlpha = 0.0;
            S.LineWidth = 1.0;
            patch(S, 'Parent', ax)
            axis(ax, 'image');  axis 'xy'
            hold(ax, 'off');

            ax = subplot(2,3,2);
            imagesc(ax,xSupportDegs,ySupportDegs, thePSFData.vLambdaWeightedPSF);
            axis(ax, 'image');  axis 'xy'
            title('Vlambda PSF');

            ax = subplot(2,3,3);
            pcolor(ax,1:stimulusPixelSamples,1:stimulusPixelSamples, theStimulusPixelImage);
            axis(ax, 'image'); axis 'xy'
            set(ax, 'XLim', [0 stimulusPixelSamples+1], 'YLim', [0 stimulusPixelSamples+1]);
            title(sprintf('stimulusPixel (%d x%d)', stimulusPixelSamples, stimulusPixelSamples));

            ax = subplot(2,3,4);
            imagesc(ax,xSupportDegs,ySupportDegs, theOpticalStimulusPSF);
            axis(ax, 'image');  axis 'xy'
            title(sprintf('Vlambda PSF x stimulusPixel (%d x%d)', stimulusPixelSamples, stimulusPixelSamples));

           
            ax = subplot(2,3,5);
            imagesc(ax,xSupportDegs,ySupportDegs, theVisualRFmap);
            axis(ax, 'image');  axis 'xy'

            S = visualRFcontourData{iRGC};
            S.FaceVertexCData = [0.85 0.85 0.85];
            S.FaceColor = 'flat';
            S.EdgeColor = [1 0 0];
            S.FaceAlpha = 0.0;
            S.LineWidth = 1.0;

            hold(ax, 'on');
            patch(S, 'Parent', ax)
            hold(ax, 'off');


            ax = subplot(2,3,6);
            hold(ax, 'on');
            S.FaceAlpha = 0.5;
            patch(S, 'Parent', ax)
            axis(ax, 'image');  axis 'xy'
            set(ax, 'XLim', [xSupportDegs(1) xSupportDegs(end)], 'YLim', [ySupportDegs(1) ySupportDegs(end)]);
            drawnow;
        end % generateFigures

    end
end


function theMidgetRGCMosaic = generateTheCenterOnlyConnectedMRGCmosaic(mosaicParams)

    % Generate mRGC mosaic along with its input coneMosaic, and connect cones to the RF centers
    theMidgetRGCMosaic = mRGCMosaic(...
            'eccentricityDegs', mosaicParams.eccDegs, ...
            'sizeDegs', mosaicParams.sizeDegs, ...
            'chromaticSpatialVarianceTradeoff', mosaicParams.chromaticSpatialVarianceTradeoff, ...
            'maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', mosaicParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs);

end


