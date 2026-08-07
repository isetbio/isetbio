%
% RGCMosaicConstructor.visualize.StrehlRatioAsAFunctionOfDefocusAndAstigmatism
%
function StrehlRatioAsAFunctionOfDefocusAndAstigmatism(...
    StrehlRatioAsAFunctionOfDefocusAndAstigmatism, ...
    theOptimalStrehlRatioDefocusAndAstigmatismDiopters, theOptimalStrehlRatio, theOptimalStrehlRatioPSF, ...
    examinedDefocusDiopters, examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, ...
    whichEye, zernikeDataBase, subjectID,  ...
    varargin)

    p = inputParser;
    p.addParameter('figureDir', [], @(x)(isempty(x)||(ischar(x))));
    p.addParameter('darkScheme', true, @islogical);
    p.addParameter('backgroundIsTransparent', false, @islogical);
    p.addParameter('pdfFileName', '', @(x)(isempty(x)||(ischar(x))));
    p.parse(varargin{:});

    figureDir = p.Results.figureDir;
    thePDFfileName = p.Results.pdfFileName;
    backgroundIsTransparent = p.Results.backgroundIsTransparent;

    
    ff = PublicationReadyPlotLib.figureComponents('2x2 standard figure', ...
        'darkScheme', p.Results.darkScheme);
    
    % Initialize figure
    hFig = figure(); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    theLUT = brewermap(1024, '*greys');


    the3DSize(1) = numel(examinedDefocusDiopters);
    the3DSize(2) = numel(examinedObliqueAstigmatismDiopters);
    the3DSize(3) = numel(examinedVerticalAstigmatismDiopters);
    StrehlRatio3Dmap = reshape(StrehlRatioAsAFunctionOfDefocusAndAstigmatism,the3DSize);
    [theMaxStrehlRatio, idx] = max(StrehlRatioAsAFunctionOfDefocusAndAstigmatism);
    [i1, i2, i3] = ind2sub(size(StrehlRatio3Dmap), idx);

    
    % Plot the slice over the defocus x oblique astigmatism plane at the
    % optimal vertical astigmatism
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, :, i3));
    ax = theAxes{1,1};
    imagesc(ax, examinedDefocusDiopters, examinedObliqueAstigmatismDiopters, StrehlRatio2Dslice');
    hold(ax, 'on');
    plot(ax, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(1), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(2), ...
        'rx', 'LineWidth', 1.5, 'MarkerSize', 16);
    xlabel(ax, 'defocus (D)');
    ylabel(ax, sprintf('oblique astigmatism (D)\n(optimal: %2.3fD)', theOptimalStrehlRatioDefocusAndAstigmatismDiopters(2)));
    axis(ax, 'square'); axis(ax, 'xy'); 
    set(ax, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(ax, theLUT);
    colorbar(ax)

    % Plot the slice over the vertical astigmatism x defocus plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, i2, :));
    ax = theAxes{2,1};
    imagesc(ax, examinedVerticalAstigmatismDiopters, examinedDefocusDiopters, StrehlRatio2Dslice);
    hold(ax, 'on');
    plot(ax, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(3), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(1), ...
        'rx', 'LineWidth', 1.5, 'MarkerSize', 16);
    xlabel(ax, 'vertical astigmatism (D)');
    ylabel(ax, sprintf('defocus (D)\n(optimal: %2.3fD)', theOptimalStrehlRatioDefocusAndAstigmatismDiopters(1)));
    axis(ax, 'square'); axis(ax, 'xy'); 
    set(ax, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(ax, theLUT);
    colorbar(ax)


    % Plot the slice over the oblique x vertical astigmatism plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(i1, :, :));
    ax = theAxes{1,2};
    imagesc(ax,examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, StrehlRatio2Dslice');
    hold(ax, 'on');
    plot(ax, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(2), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(3), ...
        'rx', 'LineWidth', 1.5, 'MarkerSize', 16);
    xlabel(ax, 'oblique astigmatism (D)');
    ylabel(ax, sprintf('vertical astigmatism (D)\n(optimal: %2.3fD)', theOptimalStrehlRatioDefocusAndAstigmatismDiopters(3)));
    axis(ax, 'square'); axis(ax, 'xy'); 
    set(ax, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(ax, theLUT);
    colorbar(ax)
    
    
    XLimsArcMin = 5*[-1 1];
    YLimsArcMin = XLimsArcMin;
    XTicksArcMin = -6:2:6;
    YTicksArcMin = -6:2:6;
 
    [~,targetInFocusWavelengthIndex] = min(abs(theOptimalStrehlRatioPSF.supportWavelength-550));

    ax = theAxes{2,2};
    RGCMosaicConstructor.helper.optics.visualizePSFatWavelength(ax, ...
        theOptimalStrehlRatioPSF, targetInFocusWavelengthIndex, ...
	    max(theOptimalStrehlRatioPSF.data(:)), ...
		sprintf('%s (Strehl ratio: %2.3f)',  whichEye, theOptimalStrehlRatio), ...
		'XLimsArcMin', XLimsArcMin, ...
		'YLimsArcMin', YLimsArcMin, ...
		'XTicksArcMin', XTicksArcMin, ... 
		'YTicksArcMin', YTicksArcMin);
    colormap(ax, theLUT);
    colorbar(ax);

    % Finalize figure using the Publication-Ready format
    ff.box = 'on';
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{1,2},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{2,1},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{2,2},ff);

    if (backgroundIsTransparent)
        set(hFig, 'Color', 'none');
        set(ax, 'Color', 'none', 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    end

    if (isempty(thePDFfileName))
        thePDFfileName = sprintf('StrehlRatio3Doptimization_%s_%s_subjID_%d', whichEye, zernikeDataBase, subjectID);
    end

    thePDFfileName = fullfile(figureDir,thePDFfileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end

