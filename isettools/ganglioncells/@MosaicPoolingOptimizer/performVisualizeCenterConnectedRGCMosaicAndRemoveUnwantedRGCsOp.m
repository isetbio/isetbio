function performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('identifyPooledCones', true, @islogical);
    p.addParameter('identifyInputCones',  true, @islogical);
    p.addParameter('plotRFoutlines',  true, @islogical);
    p.addParameter('backgroundColor', [0 0 0], @(x)(isnumeric(x)||(numel(x)==3)));
    p.parse(varargin{:});

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

    identifyPooledCones = p.Results.identifyPooledCones;
    identifyInputCones = p.Results.identifyInputCones;
    plotRFoutlines = p.Results.plotRFoutlines;
    backgroundColor = p.Results.backgroundColor;

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1]);

    theMidgetRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'identifiedConeApertureThetaSamples', 32, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones',identifyInputCones, ...
            'plotRFoutlines', plotRFoutlines, ...
            'domainVisualizationLimits', [], ...
            'domainVisualizationTicks', [], ...
            'backgroundColor', backgroundColor);

    % Report stats of how many RGCs contain each # of cones in their RFcenter
    centerConesNumCases = theMidgetRGCMosaic.centerConePoolingStats();

    % Ask user whether to eliminate certain RGCs
    eliminateRGCs = input('Eliminate RGCs with a certain no of center cones? [y = YES] : ', 's');
    if (strcmpi(eliminateRGCs, 'y'))
        % Ask user which RGCs to eliminate
        targetCenterConesNum = input(sprintf('What number of center cones [%d .. %d]? :', min(centerConesNumCases), max(centerConesNumCases)));
        
        % Eliminate the target RGCs
        theMidgetRGCMosaic.eliminateRGCsWithThisManyCenterConesNum(targetCenterConesNum);
        
        % Ask user whether to save the updated RGC mosaic
        exportUpdatedRGCMosaic = input('Export the updated RGC mosaic [y = YES] : ', 's');
        if (strcmpi(exportUpdatedRGCMosaic, 'y'))
            save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic', '-v7.3');
            fprintf('The updated RGC mosaic was saved to %s.\n', ...
                fullfile(resourcesDirectory,mosaicFileName));
        else
            fprintf('The updated RGC mosaic was not saved to the disk.\n');
        end
     else
            fprintf('No change in the RGC mosaic\n');
    end

end
