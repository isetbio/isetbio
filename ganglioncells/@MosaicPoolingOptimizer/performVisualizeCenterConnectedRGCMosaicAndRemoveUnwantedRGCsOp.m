function performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('identifyPooledCones', true, @islogical);
    p.addParameter('identifyInputCones',  true, @islogical);
    p.addParameter('plotRFoutlines',  true, @islogical);
    p.addParameter('labelRetinalMeridians', false, @islogical);
    p.addParameter('backgroundColor', [0 0 0], @(x)(isnumeric(x)||(numel(x)==3)));
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
    p.addParameter('reverseXDir', false, @islogical);
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
    labelRetinalMeridians = p.Results.labelRetinalMeridians;
    backgroundColor = p.Results.backgroundColor;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    reverseXDir = p.Results.reverseXDir;

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x1 small wide');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);

    theMidgetRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle',theAxes{1,1}, ...
            'labelRetinalMeridians', labelRetinalMeridians, ...
            'identifiedConeApertureThetaSamples', 32, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones',identifyInputCones, ...
            'plotRFoutlines', plotRFoutlines, ...
            'centerSubregionContourSamples', 32, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'backgroundColor', backgroundColor, ...
            'fontSize', ff.fontSize, ...
            'fontAngle', ff.axisFontAngle);

    if (reverseXDir)
        set(theAxes{1,1}, 'XDir', 'reverse');
    end
    
    % ============ Now generate separate figs for each parameter (PLOS format)
    rawFiguresRoot = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode/Raw';
    pdfFileName = fullfile(rawFiguresRoot, sprintf('RFcenters_%2.1fdegs.pdf', mosaicParams.eccDegs(1)));

    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    % Generate paper-ready figures (scaled versions of the figures i
    % nrawFiguresRoot directory) which are stored in the PaperReady folder
    PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
    commandString = sprintf('%s/cpdf -args %s/generatePLOSOnePaperReadyFigures.txt', PLOSdirectory, PLOSdirectory);
    system(commandString);

    pause
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
