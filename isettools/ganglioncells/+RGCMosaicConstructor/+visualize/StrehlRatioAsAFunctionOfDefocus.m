%
% RGCMosaicConstructor.visualize.StrehlRatioAsAFunctionOfDefocus
%
function StrehlRatioAsAFunctionOfDefocus(...
    examinedRefractionErrorDiopters, StrehlRatioAsAFunctionOfDefocus, ...
    theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio, ...
    whichEye, zernikeDataBase, subjectID, ...
    varargin)

    p = inputParser;
    p.addParameter('axesToRenderIn', [], @(x)(isempty(x)||(isa(x, 'handle'))));
    p.addParameter('darkScheme', false, @islogical);
    p.parse(varargin{:});
    axesToRenderIn = p.Results.axesToRenderIn;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure', ...
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
    theTitle = sprintf('%s_%s_subjID_%d.pdf', whichEye, zernikeDataBase, subjectID);
    title(ax, theTitle, 'Interpreter', 'none');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    if (isempty(axesToRenderIn))
        pdfFileName = sprintf('StrehlOptimization_%s.pdf', theTitle);
        thePDFfileName = fullfile(pdfFileName);
        NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
    end

end

