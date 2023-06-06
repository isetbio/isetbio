function generateMidgetRGCPosterFigures()

    mosaicEcc = 2.5;
    mosaicEcc = 7.0;

    % Get mosaic ecc and size
    mosaicParams = getMosaicParams(mosaicEcc);

    % Generate mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the mosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    % Detailed view shows the connections, non-detailed just the RF outline
    showDetailedNarrowView = true;

    if (showDetailedNarrowView)
        yMin =  1.0;
        yMax =  1.5;
        identifyPooledCones = true;
        identifyInputCones = true;
        %figureFormat = MSreadyPlot.figureFormat('1x1 very long poster');
        figureFormat = MSreadyPlot.figureFormat('1x1 poster');
        dpi = 1200;
        pdfFileName = sprintf('detailMosaicAt%2.1fdegs.pdf', mosaicEcc);
        pngFileName = sprintf('detailMosaicAt%2.1fdegs.png', mosaicEcc);
    else
        yMin =  -1.5;
        yMax =   1.5;
        identifyPooledCones = ~true;
        identifyInputCones = ~true;
        figureFormat = MSreadyPlot.figureFormat('1x1 poster');
        dpi = 600;
        pdfFileName = sprintf('fullMosaicAt%2.1fdegs.pdf', mosaicEcc);
        pngFileName = sprintf('fullMosaicAt%2.1fdegs.png', mosaicEcc);
    end

    xMin =  mosaicEcc-3;
    xMax =  mosaicEcc+3;
    domainVisualizationLimits = [xMin xMax yMin yMax];
    domainVisualizationTicks = struct(...
        'x', [0:0.5:10], ...
        'y', [-2:0.5:2]);

    hFig = figure(1); clf;
    MSreadyPlot.renderPosterMosaic(hFig, theMidgetRGCMosaic, '', figureFormat, ...
        'backgroundColor', [1 1 1], ...
        'identifyPooledCones', identifyPooledCones , ...
        'identifyInputCones', identifyInputCones, ...
        'domainVisualizationLimits', domainVisualizationLimits, ...
        'domainVisualizationTicks', domainVisualizationTicks);

    NicePlot.exportFigToPDF(pdfFileName, hFig, dpi);
    %NicePlot.exportFigToPNG(pngFileName, hFig, dpi);

end

function mosaicParams = getMosaicParams(mosaicEcc)

    switch (mosaicEcc)
        case 0
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
        % which covers the [1 - 4] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [0 0], ...
            'sizeDegs', [3 3]);

        case 2.5
        % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
        % which covers the [1 - 4] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [2.5 0], ...
            'sizeDegs', [3 3]);

        case 7.0
        % Mosaic params to employ. This is for the 7.0 deg - centered mosaic
        % which covers the [4-10] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [7 0], ...
            'sizeDegs', [6 3]);

        otherwise
            error('No data for this eccentricity')
    end
end

