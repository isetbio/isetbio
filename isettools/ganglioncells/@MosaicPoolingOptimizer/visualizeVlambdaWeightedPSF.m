function visualizeVlambdaWeightedPSF(theComputeReadyMRGCmosaic, opticsParams, varargin)
    
    p = inputParser;
    p.addParameter('axesHandle', [],  @(x)(isempty(x)||ishandle(x)));
    p.addParameter('figureFormat', [],  @(x)(isempty(x)||isstruct(x)));
    p.parse(varargin{:});

    axesHandle = p.Results.axesHandle;
    figureFormat = p.Results.figureFormat;

    % Generate the Vlambda weighted psfData
    thePSFData = MosaicPoolingOptimizer.generateVlambdaWeightedPSFData(...
            theComputeReadyMRGCmosaic, opticsParams);

    % Plot the vLambda weigted PSF
    tickSeparationArcMin = 3;
    psfSupportXdegs = thePSFData.supportX/60;
    psfSupportYdegs = thePSFData.supportY/60;
    psfRangeArcMin = 0.5*tickSeparationArcMin * 4;

    if (isempty(axesHandle))
        hFig = figure(1000); clf;
        figureFormat = MSreadyPlot.figureFormat('1x1 small');
        theAxes = MSreadyPlot.generateAxes(hFig,figureFormat);
        theAxes = theAxes{1,1};
    else
        theAxes = axesHandle;
        if isempty(figureFormat)
            error('No figureFormat was passed although an axesHandle was passed. If an axesHandle is passed, a figureFormat must also be passed')
        end
    end

    MSreadyPlot.render2DPSF(theAxes, ...
        psfSupportXdegs, psfSupportYdegs, ...
        thePSFData.vLambdaWeightedPSF/max(thePSFData.vLambdaWeightedPSF(:)), ...
        psfRangeArcMin/60, sprintf('V_{\\lambda}-weighted PSF'), figureFormat, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'noYLabel', ~true, ...
        'gridlessPSF', ~true, ...
        'colorMap', brewermap(1024, 'greys'));

    drawnow;
end
