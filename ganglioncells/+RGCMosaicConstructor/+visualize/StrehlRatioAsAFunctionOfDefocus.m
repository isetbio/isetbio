%
% RGCMosaicConstructor.visualize.StrehlRatioAsAFunctionOfDefocus
%
function StrehlRatioAsAFunctionOfDefocus(...
    examinedRefractionErrorDiopters, StrehlRatioAsAFunctionOfDefocus, ...
    theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio, theOptimalStrehlRatioPSF, ...
    whichEye, zernikeDataBase, subjectID,  ...
    varargin)

    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.addParameter('figureDir', [], @(x)(isempty(x)||(ischar(x))));
    p.addParameter('darkScheme', true, @islogical);
    p.addParameter('backgroundIsTransparent', false, @islogical);
    p.addParameter('pdfFileName', '', @(x)(isempty(x)||(ischar(x))));
    p.parse(varargin{:});

    axesToRenderIn = p.Results.axesToRenderIn;
    figureDir = p.Results.figureDir;
    thePDFfileName = p.Results.pdfFileName;
    backgroundIsTransparent = p.Results.backgroundIsTransparent;
    
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard very tall figure', ...
        'darkScheme', p.Results.darkScheme);
    
    if (isempty(axesToRenderIn))
        % Initialize figure
        hFig = figure(); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};
    else
        ax = axesToRenderIn;
    end




    plot(ax, examinedRefractionErrorDiopters, StrehlRatioAsAFunctionOfDefocus, 'ro-', ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0. 0.]);
    hold(ax, 'on');
    plot(ax, theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio, 'ks', ...
        'MarkerSize', 20, 'LineWidth', 3);
    plot(ax, theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio, 'cs', ...
        'MarkerSize', 20, 'LineWidth', 1.5);
  
    set(ax, 'XTick', -10:1:10, 'YLim', [0 1], 'YTick', 0:0.2:1.0);

    xlabel(ax, 'defocus (D)');
    ylabel(ax, 'Strehl ratio');
    theTitle = sprintf('%s_%s_subjID_%d\nStrehl ratio: %2.3f at %2.3fD', whichEye, zernikeDataBase, subjectID, theOptimalStrehlRatio, theOptimalStrehlRatioDefocusDiopters);
    thePDFfileName = sprintf('StrehlRatioOptimization_%s_%s_subjID_%d', whichEye, zernikeDataBase, subjectID);
    
    title(ax, theTitle, 'Interpreter', 'none');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    if (backgroundIsTransparent)
        set(hFig, 'Color', 'none');
        set(ax, 'Color', 'none', 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    end

     if (isempty(thePDFfileName))
        thePDFfileName = sprintf('StrehlRatioOptimization_%s_%s_subjID_%d', whichEye, zernikeDataBase, subjectID);
    end

    if (isempty(axesToRenderIn))
        thePDFfileName = fullfile(figureDir,thePDFfileName);
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end


    if (isempty(axesToRenderIn))
        % Also visualize the derived PSF

        % Initialize figure
        hFig = figure(); clf;
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        XLimsArcMin = 5*[-1 1];
        YLimsArcMin = XLimsArcMin;
        XTicksArcMin = -6:2:6;
        YTicksArcMin = -6:2:6;

        [~,targetInFocusWavelengthIndex] = min(abs(theOptimalStrehlRatioPSF.supportWavelength-550));
        
        RGCMosaicConstructor.helper.optics.visualizePSFatWavelength(ax, ...
            theOptimalStrehlRatioPSF, targetInFocusWavelengthIndex, ...
	        max(theOptimalStrehlRatioPSF.data(:)), ...
		    sprintf('%s (Strehl ratio: %2.3f)',  whichEye, theOptimalStrehlRatio), ...
		    'XLimsArcMin', XLimsArcMin, ...
		    'YLimsArcMin', YLimsArcMin, ...
		    'XTicksArcMin', XTicksArcMin, ... 
		    'YTicksArcMin', YTicksArcMin);
        theLUT = brewermap(1024, '*greys');
        colormap(ax, theLUT);
        colorbar(ax);

        % Finalize figure using the Publication-Ready format
        ff.box = 'on';
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        if (backgroundIsTransparent)
            set(hFig, 'Color', 'none');
            set(ax, 'Color', 'none', 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
        end
    
        if (isempty(axesToRenderIn))
            thePDFfileName = strrep(thePDFfileName, 'StrehlRatioOptimization', 'StrehlRatioOptimizationPSF');
            NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
        end

    end


end

