function retrieveAndVisualizeSubspaceRFmapForTargetRGC(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
            pdfFileName, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);

    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;

    
    if (isempty(targetCenterConeMajorityType))
            theVisualizedRGCindex = input('Enter index of target RGC : ');
            [~, theCenterConeTypeNum, targetCenterConeMajorityType, ~] = theComputeReadyMRGCmosaic.centerConeTypeWeights(theVisualizedRGCindex);
            targetCenterConesNum = sum(theCenterConeTypeNum);
            targetRGCposition = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(theVisualizedRGCindex,:);
    else
        % Find the target RGC to be visualized
        [targetCenterConesNumNotMatched, theVisualizedRGCindex] = theComputeReadyMRGCmosaic.indexOfRGCNearPosition( ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);
        if (targetCenterConesNumNotMatched)
            fprintf(2, 'Could not find an RGC with %d center cones\n', targetCenterConesNum);
            return;
        end
    end

    % Load the optimall mapped visual RF maps for all cells
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps');

    if (isempty(optimallyMappedVisualRFmaps{theVisualizedRGCindex}))
        fprintf(2, 'Optimally mapped visual RF map data for this RGC were not found in %s\n', optimallyMappedSubspaceRFmapsFileName);
        return;
    end

    % Figure format
    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x4 RF poster');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    % Plot the visual RF map
    retinalRGCRFposDegs = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(theVisualizedRGCindex,:);


    MosaicPoolingOptimizer.visualizeVisualRFmap(...
        optimallyMappedVisualRFmaps{theVisualizedRGCindex}, ...
        retinalRGCRFposDegs, ...
        theAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'reverseXDir', reverseXDir, ...
        'gridlessLineWeightingFuncions', gridlessLineWeightingFuncions, ...
        'withFigureFormat', ff);

    
    % Export to PDF
    if (isnan(targetCenterConeMajorityType))
        pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM', ...
                    targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
    else
        switch (targetCenterConeMajorityType)
           case cMosaic.LCONE_ID
                pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_LconeDominated.pdf', ...
                        targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
           case cMosaic.MCONE_ID
                pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_MconeDominated.pdf', ...
                        targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
                    
       end
    end
    pdfFileName = strrep(pdfFileName, '.pdf', pdfPostFix);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
    
end

