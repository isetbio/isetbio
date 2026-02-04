% Script to visualize full-scale position lattices
% This script is also used to generate materials for the validation
% figures for the PLOS2024 paper
%
% Usage:
%{
    t_testPositionLatticeVisualization    
%}


DegsToMMsConversionFunction = @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x);
MMsToDegsConversionFunction = @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x);


% Which mRGC/cone mosaic lattice to use
sourceLatticeSizeDegs = 64; 
whichEye = 'right eye';
neuronType = 'midget ganglion cells';  % choose between {'cones', 'midget ganglion cells'}

% Configure algorithm params
params = retinalattice.configure(sourceLatticeSizeDegs, neuronType, whichEye);
    
% Load patch generation data
load(fullfile(params.latticeGalleryDir, params.patchSaveFileName), 'dataOut', 'params', 'fovDegs', 'neuronType', 'whichEye');
      
% Unpack data
rfPositions = dataOut.rfPositions;
rfPositionsHistory = dataOut.rfPositionsHistory;
iterationsHistory = dataOut.iterationsHistory;
maxMovements = dataOut.maxMovements;
reTriangulationIterations = dataOut.reTriangulationIterations;
terminationReason = dataOut.terminationReason;
 
triangularizationsNum = size(reTriangulationIterations,2);
savedIterationsNum = size(iterationsHistory,2);
totalIterations = size(maxMovements,2);
rfsNum = size(rfPositionsHistory,2);

switch (neuronType)
    case 'cones'
        iter = savedIterationsNum;
    case 'midget ganglion cells'
        iter = savedIterationsNum-5; % 4;
end

recomputePositions = ~true;
compute2DdensityMap = ~true;

factorToMatchWatsonISETBioPeakConeDensity = 0.74;


maxRadialEccDegs = 30;
fName = sprintf('%ddegField_%s_iter_%d', 2*maxRadialEccDegs,neuronType, iter);

if (recomputePositions)
    removeOverlappingPositions = true;
    theRFpositionsMicrons = double(squeeze(rfPositionsHistory(iter,:,:)));

    % Reverse the polarity
    theRFpositionsMicrons = -theRFpositionsMicrons;

    [theRFpositionsMicrons, theRFspacingMicrons] = ...
        cropPositions(theRFpositionsMicrons, DegsToMMsConversionFunction(maxRadialEccDegs)*1e3, removeOverlappingPositions);
    theRFpositionsDegs = MMsToDegsConversionFunction(theRFpositionsMicrons * 1e-3);
    save(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons', '-v7.3');
else
    load(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons');
end

showDensityOnTopOfPositions = false;

if (~showDensityOnTopOfPositions)
    % Depict mosaic at 20 degs
    centerDegs = [-20 0];
    widthDegs = 2;
    heightDegs = 1.25;


    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', centerDegs, ...
            'width', widthDegs, ...
            'height', heightDegs, ...
            'rotation', 0.0...
        ));

    visualizedRFindices = theROI.indicesOfPointsInside(theRFpositionsDegs);
    plotMosaicPatch(100, neuronType, sprintf('%s_Ecc_%2.0fdegs', fName, centerDegs(1)), ...
        centerDegs, widthDegs, heightDegs, ...
        theRFpositionsDegs(visualizedRFindices,:));


    % Depict mosaic at 0 degs
    centerDegs = [0 0];
    widthDegs = 2;
    heightDegs = 1.25;

    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', centerDegs, ...
            'width', widthDegs, ...
            'height', heightDegs, ...
            'rotation', 0.0...
        ));

    visualizedRFindices = theROI.indicesOfPointsInside(theRFpositionsDegs);
    plotMosaicPatch(101, neuronType, sprintf('%s_Ecc_%2.0fdegs', fName, centerDegs(1)), ...
        centerDegs, widthDegs, heightDegs, ...
        theRFpositionsDegs(visualizedRFindices,:));


end



if (compute2DdensityMap)

    eccentricityDegs = [0 0]; 
    sizeDegs = [64 64];

    % sampling for the density map
    densityMapSampleSpacingDegs = 1.0;

    densityMap = computeMRGCRFmapDensityMap(theRFpositionsMicrons, theRFpositionsDegs, ...
        eccentricityDegs, sizeDegs, densityMapSampleSpacingDegs);

    mapLog = log10(densityMap.map);

    % density isocontour levels
    densityRange = prctile(mapLog(:), [1 100]);
    zLevels = 0.5:0.125:4.0;

    densityTicks = zLevels(5:8:end);
    densityTickLabels = sprintf('1E%1.0f\n',densityTicks);

    % LUT for contour plot
    theLUT = brewermap(1024, '*spectral');
    contourLineColor = [0 0 0];

    ff = PublicationReadyPlotLib.figureComponents('1x1 giant rectangular-wide mosaic');
	hFig = figure(9999); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    contourf(ax, densityMap.spatialSupportXdegs, densityMap.spatialSupportYdegs, mapLog, ...
        zLevels, 'LineColor', contourLineColor);
    axis(ax, 'xy'); axis(ax, 'image');
    set(ax, 'XLim', 32*[-1 1], 'YLim', 32*[-1 1], 'XTick', -30:5:30, 'YTick', -30:5:30);
    set(ax, 'Color', theLUT(1,:), 'CLim', densityRange, 'FontSize', 16);
    xtickangle(ax, 0);

    xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic');
    ylabel(ax, 'eccentricity, y (degs)', 'FontAngle', 'italic');
    cb = colorbar(ax, 'South', ...
        'Ticks', densityTicks, 'TickLabels', densityTickLabels);
    if (strcmp(neuronType,'cones'))
        cb.Title.String = 'cone density (cones/deg^2)';
    else
        cb.Title.String = 'mRGC RF density (RFs/deg^2)';
    end

    colormap(ax, theLUT);

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
   
    drawnow;
    NicePlot.exportFigToPDF('densityMap.pdf', hFig, 300);

end



if (~recomputePositions)
    eccDegs = sqrt(sum(theRFpositionsDegs.^2,2));
    theRFspacingDegs = MMsToDegsConversionFunction(theRFspacingMicrons* 1e-3);
    densityPerDeg2 = RGCmodels.Watson.convert.spacingToDensityForHexGrid(reshape(theRFspacingDegs, size(eccDegs)));
    
    theROIright = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [maxRadialEccDegs/4+0.15 0], ...
            'width', maxRadialEccDegs/2, ...
            'height', 0.1, ...
            'rotation', 0.0...
        ));
    theROIleft = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [-maxRadialEccDegs/4-0.15 0], ...
            'width', maxRadialEccDegs/2, ...
            'height', 0.1, ...
            'rotation', 0.0...
        ));

    theROIrightFar = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [3*maxRadialEccDegs/4 0], ...
            'width', maxRadialEccDegs, ...
            'height', 1, ...
            'rotation', 0.0...
        ));
    theROIleftFar = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [-3*maxRadialEccDegs/4 0], ...
            'width', maxRadialEccDegs, ...
            'height', 1, ...
            'rotation', 0.0...
        ));

    theROIcenter = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [0 0], ...
            'width', 0.3, ...
            'height', 0.07, ...
            'rotation', 0.0...
        ));

    visualizedRFindicesRight = theROIright.indicesOfPointsInside(theRFpositionsDegs);
    visualizedRFindicesLeft = theROIleft.indicesOfPointsInside(theRFpositionsDegs);
    visualizedRFindicesRightFar = theROIrightFar.indicesOfPointsInside(theRFpositionsDegs);
    visualizedRFindicesLeftFar = theROIleftFar.indicesOfPointsInside(theRFpositionsDegs);
    visualizedRFindicesCenter = theROIcenter.indicesOfPointsInside(theRFpositionsDegs);
    visualizedRFindices = unique([visualizedRFindicesRightFar(:); visualizedRFindicesLeftFar(:); visualizedRFindicesRight(:); visualizedRFindicesLeft(:);  visualizedRFindicesCenter(:)]);
    %visualizedRFindices = unique([visualizedRFindicesRight(:); visualizedRFindicesLeft(:)]); 
    xPos = theRFpositionsDegs(visualizedRFindices,1);
    densities = densityPerDeg2(visualizedRFindices);

    theoreticalNasalRetinaEccDegs = (xPos(xPos>=0))';
    theoreticalTemporalRetinaEccDegs = (xPos(xPos<=0))';

    rightEyeVisualFieldCorrespondingToRightNasalRetina = RGCmodels.Watson.constants.temporalMeridian;
    rightEyeVisualFieldCorrespondingToRightTemporalRetina =  RGCmodels.Watson.constants.nasalMeridian;

    switch (neuronType)
        case 'cones'
            elementSizeFactor = 1;
            [~,~,theoreticalDensityNasalRetina] = ...
                RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(...
                    theoreticalNasalRetinaEccDegs, rightEyeVisualFieldCorrespondingToRightNasalRetina, true);
            [~,~,theoreticalDensityTemporalRetina] = ...
                RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(...
                    -theoreticalTemporalRetinaEccDegs, rightEyeVisualFieldCorrespondingToRightTemporalRetina, true);
        case 'midget ganglion cells'
            elementSizeFactor = 1;
            [~,~,theoreticalDensityNasalRetinaONandOFFtypes] = ... 
                RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(...
                    theoreticalNasalRetinaEccDegs, rightEyeVisualFieldCorrespondingToRightNasalRetina);
            [~,~,theoreticalDensityTemporalRetinaONandOFFtypes] = ...
                RGCmodels.Watson.compute.midgetRGCRFSpacingAlongMeridianInRightEyeVisualField(...
                    -theoreticalTemporalRetinaEccDegs, rightEyeVisualFieldCorrespondingToRightTemporalRetina);
            % Only one subtype (ON or OFF)
            theoreticalDensityNasalRetina = 0.5*theoreticalDensityNasalRetinaONandOFFtypes;
            theoreticalDensityTemporalRetina = 0.5*theoreticalDensityTemporalRetinaONandOFFtypes;
    end

    
    if (showDensityOnTopOfPositions)
        hFig = figure(iter*10+1); clf;
        set(hFig, 'Position', [10 10 1850 1150], 'Color', [1 1 1]);
        ax = subplot('Position', [0.04 0.06 0.95 0.9]);

        scatter(ax,theRFpositionsDegs(:,1), theRFpositionsDegs(:,2), elementSizeFactor*(1.5)^2, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold(ax, 'on');
        scatter(ax,xPos, -0.95*maxRadialEccDegs*0.5 + densities*factorToMatchWatsonISETBioPeakConeDensity/17000*maxRadialEccDegs*0.9, 6^2, ...
            'MarkerFaceColor', [0 0.8 0.0], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
        scatter(ax, theoreticalNasalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityNasalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        scatter(ax, theoreticalTemporalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityTemporalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        
        set(ax, 'XLim', maxRadialEccDegs*[-1 1], 'YLim',  maxRadialEccDegs*0.5*[-1 1], 'FontSize', 20, 'FontAngle', 'italic')
        set(ax, 'XTick', -40:2:40, 'YTick', -40:2:40);
        xtickangle(ax, 0);
        axis(ax,'square'); axis(ax, 'equal'); 
        xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic')
        ylabel(ax, 'eccentricity,y (degs)', 'FontAngle', 'italic')
        drawnow;
       
        NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);
    else
        plotDensities(iter*10+1, maxRadialEccDegs, fName, neuronType, ...
            xPos, densities*factorToMatchWatsonISETBioPeakConeDensity, ...
            theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
            theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina);
    end
end


if (recomputePositions)
    removeOverlappingPositions = ~true;
    [theRFpositionsMicrons, theRFspacingMicrons] = ...
        cropPositions(theRFpositionsMicrons, DegsToMMsConversionFunction(maxRadialEccDegs)*1e3, removeOverlappingPositions);
    theRFpositionsDegs = MMsToDegsConversionFunction(theRFpositionsMicrons * 1e-3);
    save(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons', '-v7.3');
else
    load(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons');
end

if (~recomputePositions)
    
    maxRadialEccDegs = 20;
    fName = sprintf('%ddegField_%s_iter_%d', 2*maxRadialEccDegs,neuronType, iter);

    if (showDensityOnTopOfPositions)
        hFig = figure(iter*10+2); clf;
        set(hFig, 'Position', [10 10 1850 1150], 'Color', [1 1 1]);
        ax = subplot('Position', [0.04 0.06 0.95 0.9]);
        scatter(ax,theRFpositionsDegs(:,1), theRFpositionsDegs(:,2), elementSizeFactor*(2)^2, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold(ax, 'on');
        scatter(ax,xPos, -0.95*maxRadialEccDegs*0.5 + densities*factorToMatchWatsonISETBioPeakConeDensity/17000*maxRadialEccDegs*0.9, 6^2, ...
            'MarkerFaceColor', [0 0.8 0.0], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
        scatter(ax, theoreticalNasalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityNasalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        scatter(ax, theoreticalTemporalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityTemporalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        
        set(ax, 'XLim', maxRadialEccDegs*[-1 1], 'YLim',  maxRadialEccDegs*0.5*[-1 1], 'FontSize', 20, 'FontAngle', 'italic')
        set(ax, 'XTick', -40:2:40, 'YTick', -40:2:40);
        xtickangle(ax, 0);
        axis(ax,'square'); axis(ax, 'equal'); 
        xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic')
        ylabel(ax, 'eccentricity,y (degs)', 'FontAngle', 'italic')
        drawnow;
        NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);

    else

        plotDensities(iter*10+2, maxRadialEccDegs, fName, neuronType, ...
            xPos, densities*factorToMatchWatsonISETBioPeakConeDensity, ...
            theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
            theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina);
    end
end

maxRadialEccDegs = 10;
fName = sprintf('%ddegField_%s_iter_%d', 2*maxRadialEccDegs,neuronType, iter);

if (recomputePositions)
    removeOverlappingPositions = ~true;
    [theRFpositionsMicrons, theRFspacingMicrons] = ...
        cropPositions(theRFpositionsMicrons, DegsToMMsConversionFunction(maxRadialEccDegs)*1e3, removeOverlappingPositions);
    theRFpositionsDegs = MMsToDegsConversionFunction(theRFpositionsMicrons * 1e-3);
    save(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons', '-v7.3');
else
    load(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons');
end

if (~recomputePositions)
    if (showDensityOnTopOfPositions)
        
        hFig = figure(iter*10+3); clf;
        set(hFig, 'Position', [10 10 1850 1150], 'Color', [1 1 1]);
        ax = subplot('Position', [0.04 0.06 0.95 0.9]);
        scatter(ax,theRFpositionsDegs(:,1), theRFpositionsDegs(:,2), elementSizeFactor*(4)^2, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold(ax, 'on');
        scatter(ax,xPos, -0.95*maxRadialEccDegs*0.5 + densities*factorToMatchWatsonISETBioPeakConeDensity/17000*maxRadialEccDegs*0.9, 6^2, ...
            'MarkerFaceColor', [0 0.8 0.0], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
        scatter(ax, theoreticalNasalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityNasalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        scatter(ax, theoreticalTemporalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityTemporalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        
        set(ax, 'XLim', maxRadialEccDegs*[-1 1], 'YLim',  maxRadialEccDegs*0.5*[-1 1], 'FontSize', 20, 'FontAngle', 'italic')
        set(ax, 'XTick', -40:1:40, 'YTick', -40:1:40);
        xtickangle(ax, 0);
        axis(ax,'square'); axis(ax, 'equal'); 
        xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic')
        ylabel(ax, 'eccentricity,y (degs)', 'FontAngle', 'italic')
        drawnow;
        NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);
    else
        plotDensities(iter*10+3, maxRadialEccDegs, fName, neuronType, ...
            xPos, densities*factorToMatchWatsonISETBioPeakConeDensity, ...
            theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
            theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina);
    end
end

maxRadialEccDegs = 5;
fName = sprintf('%ddegField_%s_iter_%d', 2*maxRadialEccDegs,neuronType, iter);

if (recomputePositions)
    removeOverlappingPositions = ~true;
    [theRFpositionsMicrons, theRFspacingMicrons] = ...
        cropPositions(theRFpositionsMicrons, DegsToMMsConversionFunction(maxRadialEccDegs)*1e3, removeOverlappingPositions);
    theRFpositionsDegs = MMsToDegsConversionFunction(theRFpositionsMicrons * 1e-3);
    save(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons', '-v7.3');
else
    load(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons');
end

if (~recomputePositions)
    if (showDensityOnTopOfPositions)
        
        hFig = figure(iter*10+4); clf;
        set(hFig, 'Position', [10 10 1850 1150], 'Color', [1 1 1]);
        ax = subplot('Position', [0.04 0.06 0.95 0.9]);
        scatter(ax,theRFpositionsDegs(:,1), theRFpositionsDegs(:,2), elementSizeFactor*(6)^2, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold(ax, 'on');
        scatter(ax,xPos, -0.95*maxRadialEccDegs*0.5 + densities*factorToMatchWatsonISETBioPeakConeDensity/17000*maxRadialEccDegs*0.9, 6^2, ...
            'MarkerFaceColor', [0 0.8 0.0], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
        scatter(ax, theoreticalNasalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityNasalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        scatter(ax, theoreticalTemporalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityTemporalRetina/17000*maxRadialEccDegs*0.9, 4^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        
        set(ax, 'XLim', maxRadialEccDegs*[-1 1], 'YLim',  maxRadialEccDegs*0.5*[-1 1], 'FontSize', 20, 'FontAngle', 'italic')
        set(ax, 'XTick', -40:1:40, 'YTick', -40:1:40);
        xtickangle(ax, 0);
        axis(ax,'square'); axis(ax, 'equal'); 
        xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic')
        ylabel(ax, 'eccentricity,y (degs)', 'FontAngle', 'italic')
        NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);
        drawnow;
    else
         plotDensities(iter*10+4, maxRadialEccDegs, fName, neuronType, ...
            xPos, densities*factorToMatchWatsonISETBioPeakConeDensity, ...
            theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
            theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina);
    end

end

maxRadialEccDegs = 2.5;
fName = sprintf('%ddegField_%s_iter_%d', 2*maxRadialEccDegs,neuronType, iter);

if (recomputePositions)
    removeOverlappingPositions = ~true;
    [theRFpositionsMicrons, theRFspacingMicrons] = ...
        cropPositions(theRFpositionsMicrons, DegsToMMsConversionFunction(maxRadialEccDegs)*1e3, removeOverlappingPositions);
    theRFpositionsDegs = MMsToDegsConversionFunction(theRFpositionsMicrons * 1e-3);
    save(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons', '-v7.3');
else
    load(sprintf('%s.mat', fName), 'theRFpositionsDegs', 'theRFpositionsMicrons', 'theRFspacingMicrons');
end

if (~recomputePositions)
    if (showDensityOnTopOfPositions)
        
        hFig = figure(iter*10+5); clf;
        set(hFig, 'Position', [10 10 1850 1150], 'Color', [1 1 1]);
        ax = subplot('Position', [0.04 0.06 0.95 0.9]);
        scatter(ax,theRFpositionsDegs(:,1), theRFpositionsDegs(:,2), elementSizeFactor*(7)^2, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold(ax, 'on');
        scatter(ax,xPos, -0.95*maxRadialEccDegs*0.5 + densities*factorToMatchWatsonISETBioPeakConeDensity/17000*maxRadialEccDegs*0.9, 6^2, ...
            'MarkerFaceColor', [0 0.8 0.0], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
        scatter(ax, theoreticalNasalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityNasalRetina/17000*maxRadialEccDegs*0.9, 8^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        scatter(ax, theoreticalTemporalRetinaEccDegs, -0.95*maxRadialEccDegs*0.5 + theoreticalDensityTemporalRetina/17000*maxRadialEccDegs*0.9, 8^2, ...
            'MarkerFaceColor', [1 0.0 0.0], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
       
        set(ax, 'XLim', maxRadialEccDegs*[-1 1], 'YLim',  maxRadialEccDegs*0.5*[-1 1], 'FontSize', 20, 'FontAngle', 'italic')
        set(ax, 'XTick', -40:0.5:40, 'YTick', -40:0.5:40);
        axis(ax,'square'); axis(ax, 'equal'); 
        xtickangle(ax, 0);
        xlabel(ax, 'eccentricity, x (degs)', 'FontAngle', 'italic')
        ylabel(ax, 'eccentricity,y (degs)', 'FontAngle', 'italic')
        drawnow;
        NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);
    else
        plotDensities(iter*10+5, maxRadialEccDegs, fName, neuronType, ...
            xPos, densities*factorToMatchWatsonISETBioPeakConeDensity, ...
            theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
            theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina);
    end
end



function [theRFpositionsMicrons, theRFspacingsMicrons] = cropPositions(theRFpositionsMicrons, maxRadialEccMicrons, removeOverlappingPositions)
    idx = find(...
        (abs(theRFpositionsMicrons(:,1)) < maxRadialEccMicrons) & ...
        (abs(theRFpositionsMicrons(:,2)) < 0.5*maxRadialEccMicrons) );
    theRFpositionsMicrons = theRFpositionsMicrons(idx,:);
    theRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(theRFpositionsMicrons);
    if (removeOverlappingPositions)
        maxSeparationForDeclaringOverlap = 0.5;
        [idx, rfsToBeEliminated, overlapingRFindex] = cMosaic.identifyOverlappingRFs(0, 0, ...
                 theRFpositionsMicrons, theRFspacingsMicrons, maxSeparationForDeclaringOverlap);
        theRFpositionsMicrons = theRFpositionsMicrons(idx,:);
        theRFspacingsMicrons = theRFspacingsMicrons(idx);
    end
end

function plotMosaicPatch(figNo, neuronType, fName, centerDegs, ...
    widthDegs, heightDegs, theRFpositionsDegs)

    switch (neuronType)
        case 'cones'
            MarkerFaceColor = [255 120 150]/255;
            MarkerEdgeColor = 0*MarkerFaceColor*0.75;
        case 'midget ganglion cells'
            MarkerFaceColor = [0 0.9 0.0];
            MarkerEdgeColor = 0*[0 0.6 0.0];
    end

    sizeDegs = max([widthDegs heightDegs]);
    if (sizeDegs >= 20)
        deltaTick = 10;
    elseif (sizeDegs >= 10)
        deltaTick = 4;
    elseif (sizeDegs >= 5)
        deltaTick = 2;
    elseif (sizeDegs >= 2.5)
        deltaTick = 1;
    elseif (sizeDegs >= 1)
        deltaTick = 0.5;
    else
        deltaTick = 0.2;
    end

    elementSize = 0.35 * 1000/sqrt(size(theRFpositionsDegs,1));

    XTicks = -40:deltaTick:40;
    XLims = centerDegs(1) + 0.5*widthDegs*[-1 1];
    YLims = centerDegs(2) + 0.5*heightDegs*[-1 1];

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    scatter(ax,theRFpositionsDegs(:,1), ...
               theRFpositionsDegs(:,2), ...
               elementSize^2, ...
               'MarkerFaceColor', MarkerFaceColor, ...
               'MarkerEdgeColor', MarkerEdgeColor, ...
               'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1.0);
   
    axis(ax, 'equal');
    set(ax, 'XLim', XLims, 'YLim', YLims); 
    set(ax, 'XTick', XTicks, 'YTick', XTicks);

    xlabel(ax, 'eccentricity, x (degs)');
    ylabel(ax, 'eccentricity, y (degs)');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);
    drawnow;
    NicePlot.exportFigToPDF(sprintf('%sLattice.pdf', fName), hFig, 300);
        

end


function plot2DdensityMap(figNo)

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

        
    
end


function  plotDensities(figNo, maxRadialEccDegs, fName, neuronType, ...
        xPos, achievedDensities, ...
        theoreticalNasalRetinaEccDegs, theoreticalDensityNasalRetina, ...
        theoreticalTemporalRetinaEccDegs, theoreticalDensityTemporalRetina)

  
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    switch (neuronType)
        case 'cones'
            MarkerFaceColor = [255 120 150]/255;
            MarkerEdgeColor = MarkerFaceColor;
        case 'midget ganglion cells'
            MarkerFaceColor = [0.5 0.9 0.5];
            MarkerEdgeColor = MarkerFaceColor;
    end

    if (maxRadialEccDegs >= 25)
        XTicks = -40:5:40;
        MarkerFaceAlpha = 0.01;
        MarkerEdgeAlpha = 0.2;
    elseif (maxRadialEccDegs >= 20)
        XTicks = -40:4:40;
        MarkerFaceAlpha = 0.01;
        MarkerEdgeAlpha = 0.2;
    elseif (maxRadialEccDegs >= 10)
        XTicks = -40:2:40;
        MarkerFaceAlpha = 0.03;
        MarkerEdgeAlpha = 0.3;
    elseif (maxRadialEccDegs >= 5)
        XTicks = -40:1:40;
        MarkerFaceAlpha = 0.05;
        MarkerEdgeAlpha = 0.5;
    else
         XTicks = -40:0.5:40;
         MarkerFaceAlpha = 0.05;
         MarkerEdgeAlpha = 0.5;
    end

    xx = min(theoreticalTemporalRetinaEccDegs(:)):0.1:max(theoreticalTemporalRetinaEccDegs(:));
    idxTemporal = 0*xx;
    for i = 1:numel(xx)
        [~,idxTemporal(i)] = min(abs(theoreticalTemporalRetinaEccDegs-xx(i)));
    end

    xx = min(theoreticalNasalRetinaEccDegs(:)):0.1:max(theoreticalNasalRetinaEccDegs(:));
    idxNasal = 0*xx;
    for i = 1:numel(xx)
        [~,idxNasal(i)] = min(abs(theoreticalNasalRetinaEccDegs-xx(i)));
    end



    hold(ax, 'on');
    %p1 = scatter(ax, xPos, achievedDensities, 6^2, ...
    %    'MarkerFaceColor', MarkerFaceColor, 'MarkerEdgeColor', MarkerEdgeColor, ...
    %    'MarkerFaceAlpha', MarkerFaceAlpha, 'MarkerEdgeAlpha', MarkerEdgeAlpha);

    p1 = plot(ax, xPos, achievedDensities, '.', 'MarkerSize', 12, ...
        'MarkerFaceColor', MarkerFaceColor, 'MarkerEdgeColor', MarkerEdgeColor);

    p2 = plot(ax, theoreticalNasalRetinaEccDegs(idxNasal), theoreticalDensityNasalRetina(idxNasal), 'k-', 'LineWidth', 1.5);
    %plot(ax, theoreticalNasalRetinaEccDegs(idxNasal), theoreticalDensityNasalRetina(idxNasal), 'w--', 'LineWidth', 1.5);
    plot(ax, theoreticalTemporalRetinaEccDegs(idxTemporal), theoreticalDensityTemporalRetina(idxTemporal), 'k-', 'LineWidth', 1.5);
    %plot(ax, theoreticalTemporalRetinaEccDegs(idxTemporal), theoreticalDensityTemporalRetina(idxTemporal), 'w--', 'LineWidth', 1.5);
    XLims = maxRadialEccDegs*[-1 1]; 
    YLims = [0 17500];

    grid(ax, 'on');
    set(ax, 'XLim', XLims); 
    set(ax, 'XTick', XTicks, 'YTick', [0 2.5e3 5e3 7.5e3 10e3 12.5e3 15e3 17.5e3 20e3] , 'YLim', YLims);
    set(ax, 'YTickLabel', {'0', '', '5k', '', '10k', '', '15k', '', '20k'});
    xlabel(ax, 'eccentricity, x (degs)');
    legend(ax, [p1 p2], {'achieved', 'target'});
    switch (neuronType)
        case 'cones'
            ylabel(ax, 'density (cones/deg^2)');
        case 'midget ganglion cells'
            ylabel(ax, 'density (mRGCs/deg^2)');
    end

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);
    PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);
    drawnow;
    NicePlot.exportFigToPDF(sprintf('%s.pdf', fName), hFig, 300);
    
end




function densityMap = computeMRGCRFmapDensityMap(mRGCRFposMicrons, mRGCRFposDegs, eccentricityDegs, sizeDegs, deltaDegs)

    mRGCRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(mRGCRFposDegs);

    xPosDegs = eccentricityDegs(1):deltaDegs:(eccentricityDegs(1)+0.5*sizeDegs(1));
    xPosDegs = unique([-fliplr(xPosDegs) xPosDegs]);
    yPosDegs = eccentricityDegs(2):deltaDegs:(eccentricityDegs(2)+0.5*sizeDegs(2));
    yPosDegs = unique([-fliplr(yPosDegs) yPosDegs]);

    sampledPositionsDegs{1} = xPosDegs;
    sampledPositionsDegs{2} = yPosDegs ;    
    densityMap.map = cMosaic.densityMap(mRGCRFposDegs, mRGCRFspacingsDegs, sampledPositionsDegs);
    densityMap.spatialSupportXdegs = xPosDegs ;
    densityMap.spatialSupportYdegs = yPosDegs ;
end
