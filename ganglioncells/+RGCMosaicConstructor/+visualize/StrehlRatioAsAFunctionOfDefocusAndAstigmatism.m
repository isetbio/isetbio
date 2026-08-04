%
% RGCMosaicConstructor.visualize.StrehlRatioAsAFunctionOfDefocusAndAstigmatism
%
function StrehlRatioAsAFunctionOfDefocusAndAstigmatism(...
    examinedRefractionErrorDioptersObliqueAndVerticalAstigmatismErrorsMicronsGrid, ...
    StrehlRatioAsAFunctionOfDefocusAndAstigmatism, ...
    theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons, theOptimalStrehlRatio, ...
    examinedRefractionErrorDiopters, examinedObliqueAstigmatismErrorsMicrons, examinedVerticalAstigmatismErrorsMicrons, ...
    whichEye, zernikeDataBase, subjectID,  ...
    varargin)

    p = inputParser;
    p.addParameter('figureDir', [], @(x)(isempty(x)||(ischar(x))));
    p.addParameter('darkScheme', false, @islogical);
    p.addParameter('backgroundIsTransparent', false, @islogical);
    p.parse(varargin{:});

    figureDir = p.Results.figureDir;
    backgroundIsTransparent = p.Results.backgroundIsTransparent;

    
    StrehlRatio3Dmap = reshape(StrehlRatioAsAFunctionOfDefocusAndAstigmatism, [...
            numel(examinedRefractionErrorDiopters) ...
            numel(examinedObliqueAstigmatismErrorsMicrons) ...
            numel(examinedVerticalAstigmatismErrorsMicrons)]);

    idx = find(...
        (examinedRefractionErrorDioptersObliqueAndVerticalAstigmatismErrorsMicronsGrid(:,1) == theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(1)) & ...
        (examinedRefractionErrorDioptersObliqueAndVerticalAstigmatismErrorsMicronsGrid(:,2) == theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(2)) & ...
        (examinedRefractionErrorDioptersObliqueAndVerticalAstigmatismErrorsMicronsGrid(:,3) == theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(3)) ...
    );

    [i1, i2, i3] = ind2sub(size(StrehlRatio3Dmap), idx);


    ff = PublicationReadyPlotLib.figureComponents('2x2 standard figure', ...
        'darkScheme', p.Results.darkScheme);
    
    % Initialize figure
    hFig = figure(); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    theLUT = brewermap(1024, '*greys');

    % Plot the slice over the defocus x oblique astigmatism plane at the
    % optimal vertical astigmatism
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, :, i3));
    imagesc(theAxes{1,1}, examinedRefractionErrorDiopters, examinedObliqueAstigmatismErrorsMicrons, StrehlRatio2Dslice');
    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1}, theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(1), theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(2), 'bs', 'MarkerSize', 20);
    xlabel(theAxes{1,1}, 'defocus (D)');
    ylabel(theAxes{1,1}, 'oblique astigmatism (microns)');
    axis(theAxes{1,1}, 'square'); axis(theAxes{1,1}, 'xy'); 
    set(theAxes{1,1}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{1,1}, theLUT);
    colorbar(theAxes{1,1})

    % Plot the slice over the defocus x vertical astigmatism plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(:, i2, :));
    imagesc(theAxes{1,2},examinedRefractionErrorDiopters, examinedVerticalAstigmatismErrorsMicrons, StrehlRatio2Dslice');
    hold(theAxes{1,2}, 'on');
    plot(theAxes{1,2}, theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(1), theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(3), 'bs', 'MarkerSize', 20);
    xlabel(theAxes{1,2}, 'defocus (D)');
    ylabel(theAxes{1,2}, 'vertical astigmatism (microns)');
    axis(theAxes{1,2}, 'square'); axis(theAxes{1,2}, 'xy'); 
    set(theAxes{1,2}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{1,2}, theLUT);
    colorbar(theAxes{1,2})


    % Plot the slice over the oblique x vertical astigmatism plane
    StrehlRatio2Dslice = squeeze(StrehlRatio3Dmap(i1, :, :));
    imagesc(theAxes{2,1},examinedObliqueAstigmatismErrorsMicrons, examinedVerticalAstigmatismErrorsMicrons, StrehlRatio2Dslice');
    hold(theAxes{2,1}, 'on');
    plot(theAxes{2,1}, theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(1), theOptimalStrehlRatioDefocusDioptersObliqueAndVerticalAstigmatismErrorsMicrons(3), 'bs', 'MarkerSize', 20);
    xlabel(theAxes{2,1}, 'oblique astigmatism (microns)');
    ylabel(theAxes{2,1}, 'vertical astigmatism (microns)');
    axis(theAxes{2,1}, 'square'); axis(theAxes{2,1}, 'xy'); 
    set(theAxes{2,1}, 'CLim', [0 theOptimalStrehlRatio]);
    colormap(theAxes{2,1}, theLUT);
    colorbar(theAxes{2,1})
    
    
    title(ax, theTitle, 'Interpreter', 'none');

    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{1,2},ff);
    PublicationReadyPlotLib.applyFormat(theAxes{2,1},ff);

    if (backgroundIsTransparent)
        set(hFig, 'Color', 'none');
        set(ax, 'Color', 'none', 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    end

    thePDFfileName = sprintf('StrehlRationOptimization_%s_%s_subjID_%d', whichEye, zernikeDataBase, subjectID);
    thePDFfileName = fullfile(figureDir,thePDFfileName);
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);

end

