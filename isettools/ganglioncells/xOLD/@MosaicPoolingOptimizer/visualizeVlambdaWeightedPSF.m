function visualizeVlambdaWeightedPSF(theComputeReadyMRGCmosaic, opticsParams, varargin)
    
    p = inputParser;
    p.addParameter('axesHandle', [],  @(x)(isempty(x)||ishandle(x)));
    p.addParameter('figureFormat', [],  @(x)(isempty(x)||isstruct(x)));
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('gridlessPSF', false, @islogical);
    p.parse(varargin{:});

    axesHandle = p.Results.axesHandle;
    figureFormat = p.Results.figureFormat;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    gridlessPSF = p.Results.gridlessPSF;

    % Generate the Vlambda weighted psfData
    thePSFData = MosaicPoolingOptimizer.generateVlambdaWeightedPSFData(...
            theComputeReadyMRGCmosaic, opticsParams);

    % Plot the vLambda weigted PSF
    psfSupportXdegs = thePSFData.supportX/60;
    psfSupportYdegs = thePSFData.supportY/60;
    psfRangeArcMin = tickSeparationArcMin * 4;

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

    peakPSFamplitude = max(thePSFData.vLambdaWeightedPSF(:));
    MSreadyPlot.render2DPSF(theAxes, ...
        psfSupportXdegs, psfSupportYdegs, ...
        thePSFData.vLambdaWeightedPSF/peakPSFamplitude, ...
        psfRangeArcMin/60, sprintf('V_{\\lambda}-weighted PSF (max:%g)', peakPSFamplitude), figureFormat, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'noYLabel', true, ...
        'gridlessPSF', gridlessPSF, ...
        'colorMap', brewermap(1024, 'greys'));

    drawnow;
end
