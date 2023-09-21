function visualizeConePoolingRFmapAndVisualSTFforTargetRGC(...
    computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename, pdfFileName, ...
    targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;

    load(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic');
    load(mRGCMosaicSTFresponsesFilename, 'spatialFrequenciesTested', 'theMRGCMosaicOptimalSTFs');

    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x4 RF poster');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);

    % Plot the retinal cone pooling RF map
    theVisualizedRGCindex = theComputeReadyMRGCmosaic.visualizeRetinalConePoolingRFmapNearPosition(...
        targetRGCposition, targetCenterConesNum, ...
        targetCenterConeMajorityType, ...
        'theAxes', theAxes, ...
        'tickSeparationArcMin', tickSeparationArcMin, ...
        'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
        'withFigureFormat', ff, ...
        'reverseXDir', reverseXDir, ...
        'gridlessLineWeightingFuncions', gridlessLineWeightingFuncions);

    theVisualizedRGCvisualSTFdata = theMRGCMosaicOptimalSTFs{theVisualizedRGCindex};
    idx = find(strcmp(theVisualizedRGCvisualSTFdata.DoGfitParams.names, 'kS/kC'));
    KsKcRatio = theVisualizedRGCvisualSTFdata.DoGfitParams.finalValues(idx);
    idx = find(strcmp(theVisualizedRGCvisualSTFdata.DoGfitParams.names, 'RsToRc'));
    RsRcRatio = theVisualizedRGCvisualSTFdata.DoGfitParams.finalValues(idx);
    SCintSensRatio = KsKcRatio * (RsRcRatio)^2;

    % Now plot the visual STF
    plotTitle = sprintf('R_{srnd}/R_{cntr} = %2.2f, S_{srnd}/S_{cntr} = %2.2f', RsRcRatio, SCintSensRatio);
    MSreadyPlot.renderSTF(theAxes{1,4}, ...
         spatialFrequenciesTested, theVisualizedRGCvisualSTFdata.measured, ...
         theVisualizedRGCvisualSTFdata.DoGfit.sfHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.compositeSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.centerSTFHiRes, ...
         theVisualizedRGCvisualSTFdata.DoGfit.surroundSTFHiRes, ...
         plotTitle, [], ff, ...
         'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
         'noYLabel', ~true, ...
         'noYTickLabel', ~true);

    % Export to PDF
    if (isnan(targetCenterConeMajorityType))
        pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM', ...
                    targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
    else
        if (isempty(targetCenterConeMajorityType))
            [~, theCenterConeTypeNum, targetCenterConeMajorityType, ~] = theComputeReadyMRGCmosaic.centerConeTypeWeights(theVisualizedRGCindex);
            targetCenterConesNum = sum(theCenterConeTypeNum);
        end
        if (isempty(targetRGCposition))
            targetRGCposition = theComputeReadyMRGCmosaic.rgcRFpositionsDegs(theVisualizedRGCindex,:);
        end

        if (isnan(targetCenterConeMajorityType))
            pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);
        else
            switch (targetCenterConeMajorityType)
               case cMosaic.LCONE_ID
                    pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_LconeDominated_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);
               case cMosaic.MCONE_ID
                    pdfPostFix = sprintf('_atPosition_%2.2f_%2.2f_CenterConesNum_%d_MconeDominated_RGCindex_%d.pdf', ...
                            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, theVisualizedRGCindex);                    
            end
        end
    end

    pdfFileName = strrep(pdfFileName, '.pdf', pdfPostFix);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300); 
end