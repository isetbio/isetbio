%
%RGCMosaicAnalyzer.compute.MRGCtemporalFiltersFromPhotocurrentsBasedTTF
%
%
function MRGCtemporalFiltersFromPhotocurrentsBasedTTF(...
    stimulusShape, ...
    theTargetRGCwithIndex, ...
    theMRGCMosaicTTFResponsesFullFileName, ...
    theAnalyzedTTFsFullFileName)



    load(theMRGCMosaicTTFResponsesFullFileName, ...
        'theMRGCMosaic', 'stimParams', 'TTFparamsStruct', ...
        'computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex', ...
        'theMRGCMosaicTTFresponsesAllConditions', ...
        'theMRGCMosaicTemporalSupportSecondsAllConditions');

    temporalFrequenciesExamined = TTFparamsStruct.tfSupport

    responseRange = max(max(abs(squeeze(theMRGCMosaicTTFresponsesAllConditions(:,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex)))))*[-1 1];
    for iTF = 1:numel(temporalFrequenciesExamined)

        % Retrieve the photocurrent response data for this TF
        theTemporalSupportSecondsForThisTF = squeeze(theMRGCMosaicTemporalSupportSecondsAllConditions(iTF,:));
        theMRGCMosaicPhotocurrentResponsesForThisTF = squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,:));
        nanIndices = find(isnan(theTemporalSupportSecondsForThisTF));
        if (isempty(nanIndices))
            theTimeBins = 1:numel(theTemporalSupportSecondsForThisTF);
        else
            theTimeBins = 1:(nanIndices-1);
        end
        theTemporalSupportSecondsForThisTF = theTemporalSupportSecondsForThisTF(theTimeBins);
        theMRGCMosaicPhotocurrentResponsesForThisTF = theMRGCMosaicPhotocurrentResponsesForThisTF(theTimeBins,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex);

        hFig = figure(1);
        set(hFig, 'Position', [10 10 500 300]);
        plot(theTemporalSupportSecondsForThisTF*1e3, theMRGCMosaicPhotocurrentResponsesForThisTF, 'ko-');
        set(gca, 'XLim', [0 2000], 'YLim', responseRange);
        xlabel('time (msec)');
        ylabel('pAmps');
        title(gca, sprintf('%s, TF = %2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF)));
        drawnow;
    end % iTF
end


