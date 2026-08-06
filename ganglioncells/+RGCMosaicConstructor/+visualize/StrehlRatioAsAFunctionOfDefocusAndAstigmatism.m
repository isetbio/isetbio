%
% RGCMosaicConstructor.visualize.StrehlRatioAsAFunctionOfDefocusAndAstigmatism
%
function StrehlRatioAsAFunctionOfDefocusAndAstigmatism(...
    StrehlRatioAsAFunctionOfDefocusAndAstigmatism, ...
    theOptimalStrehlRatioDefocusAndAstigmatismDiopters, theOptimalStrehlRatio, theOptimalStrehlRatioPSF, ...
    examinedRefractionErrorDiopters, examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, ...
    whichEye, zernikeDataBase, subjectID,  ...
    varargin)

    p = inputParser;
    p.addParameter('figureDir', [], @(x)(isempty(x)||(ischar(x))));
    p.addParameter('darkScheme', true, @islogical);
    p.addParameter('backgroundIsTransparent', false, @islogical);
    p.parse(varargin{:});

    figureDir = p.Results.figureDir;
    backgroundIsTransparent = p.Results.backgroundIsTransparent;

    
    ff = PublicationReadyPlotLib.figureComponents('2x2 standard figure', ...
        'darkScheme', p.Results.darkScheme);
    
    % Initialize figure
    hFig = figure(); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    theLUT = brewermap(1024, '*greys');


    the3DSize(1) = numel(examinedRefractionErrorDiopters);
    the3DSize(2) = numel(examinedObliqueAstigmatismDiopters);
    the3DSize(3) = numel(examinedVerticalAstigmatismDiopters);
    StrehlRatio3Dmap = reshape(StrehlRatioAsAFunctionOfDefocusAndAstigmatism,the3DSize);
    [theMaxStrehlRatio, idx] = max(StrehlRatioAsAFunctionOfDefocusAndAstigmatism);
    [i1, i2, i3] = ind2sub(size(StrehlRatio3Dmap), idx);

    
    % Plot the slice over the defocus x oblique astigmatism plane at the
    % optimal vertical astigmatism
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, :, i3));
    imagesc(theAxes{1,1}, examinedRefractionErrorDiopters, examinedObliqueAstigmatismDiopters, StrehlRatio2Dslice');
    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1}, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(1), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(2), ...
        'rx', 'MarkerSize', 20);
    xlabel(theAxes{1,1}, 'defocus (D)');
    ylabel(theAxes{1,1}, 'oblique astigmatism (D)');
    axis(theAxes{1,1}, 'square'); axis(theAxes{1,1}, 'xy'); 
    set(theAxes{1,1}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{1,1}, theLUT);
    colorbar(theAxes{1,1})

    % Plot the slice over the defocus x vertical astigmatism plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, i2, :));
    imagesc(theAxes{1,2},examinedRefractionErrorDiopters, examinedVerticalAstigmatismDiopters, StrehlRatio2Dslice');
    hold(theAxes{1,2}, 'on');
    plot(theAxes{1,2}, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(1), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(3), ...
        'rx', 'MarkerSize', 20);
    xlabel(theAxes{1,2}, 'defocus (D)');
    ylabel(theAxes{1,2}, 'vertical astigmatism (D)');
    axis(theAxes{1,2}, 'square'); axis(theAxes{1,2}, 'xy'); 
    set(theAxes{1,2}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{1,2}, theLUT);
    colorbar(theAxes{1,2})


    % Plot the slice over the oblique x vertical astigmatism plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(i1, :, :));
    imagesc(theAxes{2,1},examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, StrehlRatio2Dslice');
    hold(theAxes{2,1}, 'on');
    plot(theAxes{2,1}, ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(2), ...
        theOptimalStrehlRatioDefocusAndAstigmatismDiopters(3), ...
        'rx', 'MarkerSize', 20);
    xlabel(theAxes{2,1}, 'oblique astigmatism (D)');
    ylabel(theAxes{2,1}, 'vertical astigmatism (D)');
    axis(theAxes{2,1}, 'square'); axis(theAxes{2,1}, 'xy'); 
    set(theAxes{2,1}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{2,1}, theLUT);
    colorbar(theAxes{2,1})
    
    
    XLimsArcMin = 5*[-1 1];
    YLimsArcMin = XLimsArcMin;
    XTicksArcMin = -6:2:6;
    YTicksArcMin = -6:2:6;
 
    [~,targetInFocusWavelengthIndex] = min(abs(theOptimalStrehlRatioPSF.supportWavelength-550));

    RGCMosaicConstructor.helper.optics.visualizePSFatWavelength(theAxes{2,2}, ...
        theOptimalStrehlRatioPSF, targetInFocusWavelengthIndex, ...
	    max(theOptimalStrehlRatioPSF.data(:)), ...
		sprintf('optimal Strehl ratio (%2.3f)', theOptimalStrehlRatio), ...
		'XLimsArcMin', XLimsArcMin, ...
		'YLimsArcMin', YLimsArcMin, ...
		'XTicksArcMin', XTicksArcMin, ... 
		'YTicksArcMin', YTicksArcMin);
    colormap(theAxes{2,2}, brewermap(1024, '*greys'));



    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{1,2},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{2,1},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{2,2},ff);

    if (backgroundIsTransparent)
        set(hFig, 'Color', 'none');
        set(ax, 'Color', 'none', 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    end

    thePDFfileName = sprintf('StrehlRationOptimization_%s_%s_subjID_%d', whichEye, zernikeDataBase, subjectID);
    thePDFfileName = fullfile(figureDir,thePDFfileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);

end

