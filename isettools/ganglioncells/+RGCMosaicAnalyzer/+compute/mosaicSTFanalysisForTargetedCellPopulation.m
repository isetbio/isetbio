%
%
% RGCMosaicAnalyzer.compute.mosaicSTFanalysisForTargetedCellPopulation()
%
%
function mosaicSTFanalysisForTargetedCellPopulation(...
    theMRGCMosaicSTFResponsesFullFileName, ...
    theAnalyzedSTFsFileName, ...
    targetedCenterPurityRange, ...
    targetedCenterConeNumerosityRange, ...
    targetedSurroundPurityRange, ...
    targetedRadialEccentricityRange, ...
    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses, ...
    visualizeConeExcitationVsPhotocurrentSTFs, varargin)

    p = inputParser;
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    % Execute the parser
    p.parse(varargin{:});
    exportPDFdirectory = p.Results.exportPDFdirectory;

    theData = who('-file', theMRGCMosaicSTFResponsesFullFileName);
    photocurrentBasedSTFsComputed = ismember('theMRGCMosaicPcurrentBasedResponseTemporalSupportSeconds', theData);

    theNoiseFreeSpatioTemporalMRGCMosaicPcurrentBasedResponses2DSTF = [];
    theMRGCMosaicPcurrentBasedResponseTemporalSupportSeconds = [];

    if (photocurrentBasedSTFsComputed)
        load(theMRGCMosaicSTFResponsesFullFileName, ...
            'theMRGCMosaic', 'stimParams', ...
            'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF', ...
            'theMRGCMosaicResponseTemporalSupportSeconds', ...
            'theNoiseFreeSpatioTemporalMRGCMosaicPcurrentBasedResponses2DSTF', ...
            'theMRGCMosaicPcurrentBasedResponseTemporalSupportSeconds');
    else
        load(theMRGCMosaicSTFResponsesFullFileName, ...
            'theMRGCMosaic', 'stimParams', ...
            'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF', ...
            'theMRGCMosaicResponseTemporalSupportSeconds');
    end

    allowNonZeroMeanResponse = ~true;

    [theConeModulationsBasedBPIs, thePhotocurrentsBasedBPIs, ...
     theConeModulationsBasedTargetSTFamplitudeSpectra , ...
     thePhotocurrentsBasedTargetSTFamplitudeSpectra, ...
     targetRGCindices, theCenterConeDominances, ...
     theCenterConeNumerosities, theCenterConePurities, ...
     theSurroundConePurities] = analyzeSTFresponsesForTargetMRGCs(...
        stimParams, ...
        theMRGCMosaic, ...
        theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF, ...
        theMRGCMosaicResponseTemporalSupportSeconds, ...
        theNoiseFreeSpatioTemporalMRGCMosaicPcurrentBasedResponses2DSTF, ...
        theMRGCMosaicPcurrentBasedResponseTemporalSupportSeconds, ...
        allowNonZeroMeanResponse, ...
        targetedCenterPurityRange, ...
        targetedCenterConeNumerosityRange, ...
        targetedSurroundPurityRange, ...
        targetedRadialEccentricityRange, ...
        theAnalyzedSTFsFileName, ...
        visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses, ...
        visualizeConeExcitationVsPhotocurrentSTFs, ...
        exportPDFdirectory);


    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1150 500], 'Color', [1 1 1])

    if (all(stimParams.coneContrasts == [1 0 0]))
        stimulusChroma = 'L-isolating';
    elseif (all(stimParams.coneContrasts == [0 1 0]))
        stimulusChroma = 'M-isolating';
    elseif (all(stimParams.coneContrasts == [1 1 1]))
        stimulusChroma = 'achromatic';
    else
        error('Unknown stimulus chromaticity: %2.1f %2.1f %2.1f', ...
                stimParams.coneContrasts(1), stimParams.coneContrasts(2), stimParams.coneContrasts(3));
    end

    LconeCenterIdx = find(theCenterConeDominances == cMosaic.LCONE_ID);
    MconeCenterIdx = find(theCenterConeDominances == cMosaic.MCONE_ID);

    for iOri = 1:numel(stimParams.orientationDegs)
        ax = subplot(1,numel(stimParams.orientationDegs),iOri);
        hold(ax, 'on')

        thePlottedLconeCenterIndices = [];
        thePlottedMconeCenterIndices = [];

        if (ismember(stimulusChroma, {'L-isolating', 'achromatic'}))
            if (numel(LconeCenterIdx)>0)
                idx = find(~isnan(squeeze(theConeModulationsBasedBPIs(LconeCenterIdx, iOri))));
                thePlottedLconeCenterIndices = LconeCenterIdx(idx);
                p1 = scatter(ax, ...
                     theConeModulationsBasedBPIs(thePlottedLconeCenterIndices, iOri), ...
                     thePhotocurrentsBasedBPIs(thePlottedLconeCenterIndices, iOri), 144, ...
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
                    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
            end
        end


        if (ismember(stimulusChroma, {'M-isolating', 'achromatic'}))
            if (numel(MconeCenterIdx)>0)
                idx = find(~isnan(squeeze(theConeModulationsBasedBPIs(MconeCenterIdx, iOri))));
                thePlottedMconeCenterIndices  = MconeCenterIdx(idx);
                p2 = scatter(ax, ...
                     theConeModulationsBasedBPIs(thePlottedMconeCenterIndices, iOri), ...
                     thePhotocurrentsBasedBPIs(thePlottedMconeCenterIndices, iOri), 144, ...
                    'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
                    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
            end
        end

        % Identity line
        plot(ax, [0 1], [0 1], 'k-', 'LineWidth', 1.0);

        if ((numel(thePlottedLconeCenterIndices)>0) && (numel(thePlottedMconeCenterIndices)>0))
            legend(ax, [p1 p2], {'L-cone dominated', 'M-cone dominated'}, 'Location', 'SouthEast');
        end

        if ((numel(thePlottedMconeCenterIndices)==0) && (numel(thePlottedLconeCenterIndices)>0))
            legend(ax, p1, {'L-cone dominated'}, 'Location', 'SouthEast');
        end

        if ((numel(thePlottedLconeCenterIndices)==0) && (numel(thePlottedMconeCenterIndices)>0))
            legend(ax, p2, {'M-cone dominated'}, 'Location', 'SouthEast');
        end

        xtickangle(ax, 0);
        axis(ax, 'square')
        set(ax, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1)
        grid(ax, 'on')
        title(ax, sprintf('%s (%2.0f%%)\norientation = %d degs', ...
            stimulusChroma, 100*stimParams.contrast, stimParams.orientationDegs(iOri)));
        set(ax, 'FontSize', 20);
        xlabel(ax, 'BPI (cone modulations-based STF)');
        if (iOri == 1)
            ylabel(ax, 'BPI (photocurrent-based STF)');
        end
    end % for iOri
end

function [theConeModulationsBasedBPIs, thePhotocurrentsBasedBPIs, ...
          theConeModulationsBasedTargetSTFamplitudeSpectra , ...
          thePhotocurrentsBasedTargetSTFamplitudeSpectra, ...
          targetRGCindices, theCenterConeDominances, ...
          theCenterConeNumerosities, theCenterConePurities, theSurroundConePurities] = analyzeSTFresponsesForTargetMRGCs(...
    stimParams, ...
    theMRGCMosaic, ...
    theMRGCMosaicConeModulationsBasedResponses, ...
    theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
    theMRGCMosaicPhotocurrentBasedResponses, ...
    theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
    allowNonZeroMeanResponse, ...
    targetedCenterPurityRange, ...
    targetedCenterConeNumerosityRange, ...
    targetedSurroundPurityRange, ...
    targetedRadialEccentricityRange, ...
    analyzedSTFsFileName, ...
    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses, ...
    visualizeConeExcitationVsPhotocurrentSTFs, ...
    exportPDFdirectory)

    % Find the RGCs with the desired target properties
    [targetRGCindices, theSurroundConePurities, theCenterConeDominances, ...
        theCenterConeNumerosities, theCenterConePurities] = theMRGCMosaic.indicesOfRGCsWithinTargetedPropertyRanges( ...
        targetedCenterConeNumerosityRange, ...
        targetedSurroundPurityRange, ...
        targetedRadialEccentricityRange, ...
        targetedCenterPurityRange);

    fprintf('\n\nStats for %d mRGCs (out of a total of %d)\n', numel(targetRGCindices), theMRGCMosaic.rgcsNum);
    examinedCenterConePurityRange     = [min(theCenterConePurities(:)) max(theCenterConePurities(:))];
    examinedCenterConeDominanceRange  = [min(theCenterConeDominances(:)) max(theCenterConeDominances(:))];
    examinedCenterConeNumerosityRange = [min(theCenterConeNumerosities(:)) max(theCenterConeNumerosities(:))];
    examinedSurroundConePurityRange   = [min(theSurroundConePurities(:)) max(theSurroundConePurities(:))];
    fprintf('\n\nStats for %d mRGCs (out of a total of %d)\n', numel(targetRGCindices), theMRGCMosaic.rgcsNum);
    fprintf('\n centerConePurityRange: %2.2f - %2.2f\n', examinedCenterConePurityRange(1), examinedCenterConePurityRange(2));
    fprintf('\n centerConePurityRange: %2.0f - %2.0f\n', examinedCenterConeNumerosityRange(1), examinedCenterConeNumerosityRange(2));
    fprintf('\n centerConeDominanceRange: %2.2f - %2.2f\n', examinedCenterConeDominanceRange(1), examinedCenterConeDominanceRange(2));
    fprintf('\n surroundConePurityRange: %2.2f - %2.2f\n', examinedSurroundConePurityRange(1), examinedSurroundConePurityRange(2));
    fprintf('\n\n');

    % Sort cells according to their center cone numerosity
    [~,idx] = sort(theCenterConeNumerosities, 'ascend');

    targetRGCindices = targetRGCindices(idx);
    theCenterConePurities = theCenterConePurities(idx);
    theCenterConeNumerosities = theCenterConeNumerosities(idx);
    theCenterConeDominances = theCenterConeDominances(idx);
    theSurroundConePurities = theSurroundConePurities(idx);

    orientationsNum = size(theMRGCMosaicConeModulationsBasedResponses,1);
    sfsNum = size(theMRGCMosaicConeModulationsBasedResponses,2);

    theConeModulationsBasedTargetSTFamplitudeSpectra = zeros(numel(targetRGCindices), orientationsNum, sfsNum);
    thePhotocurrentsBasedTargetSTFamplitudeSpectra = zeros(numel(targetRGCindices), orientationsNum, sfsNum);
    theConeModulationsBasedTargetSTFphaseSpectra = zeros(numel(targetRGCindices), orientationsNum, sfsNum);
    thePhotocurrentsBasedTargetSTFphaseSpectra = zeros(numel(targetRGCindices), orientationsNum, sfsNum);

    fprintf('Computing STFs for %d target mRGCs', numel(targetRGCindices));

    % The computed BPIs
    theConeModulationsBasedBPIs = nan(numel(targetRGCindices), orientationsNum);
    thePhotocurrentsBasedBPIs = nan(numel(targetRGCindices), orientationsNum);

    if (visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses)
        % Serial for loop if we are visualizing
        for iRGC = 1:numel(targetRGCindices)
            theRGCindex = targetRGCindices(iRGC);

            for iOri = 1:orientationsNum

                theConeModulationResponses = squeeze(theMRGCMosaicConeModulationsBasedResponses(iOri, :, :, theRGCindex));

                [theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 theConeModulationsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFspectra(...
                    theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                    theConeModulationResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    allowNonZeroMeanResponse, ...
                    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses);

                pause
                % Remove mean of the pCurrent responses before fitting the sinusoid
                thePCurrentResponses = squeeze(theMRGCMosaicPhotocurrentBasedResponses(iOri, :, :, theRGCindex));
                theMeanPhotoCurrentResponsesAcrossEachTimeCourse = mean(thePCurrentResponses,2);
                thePCurrentResponses = bsxfun(@minus, thePCurrentResponses, theMeanPhotoCurrentResponsesAcrossEachTimeCourse);

                [thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFspectra(...
                    theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                    thePCurrentResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    allowNonZeroMeanResponse, ...
                    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses);

               
                pause
                coneModulationsSTF = squeeze(theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));
                photocurrentsSTF = squeeze(thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));

                theConeModulationsBasedBPIs(iRGC, iOri) = ...
                    bandPassIndex(stimParams.spatialFrequencyCPD, coneModulationsSTF);

                thePhotocurrentsBasedBPIs(iRGC, iOri) = ...
                    bandPassIndex(stimParams.spatialFrequencyCPD, photocurrentsSTF);

            end % iOri
        end % for iRGC
    else
        % Do a parfor
        parfor iRGC = 1:numel(targetRGCindices)
            theRGCindex = targetRGCindices(iRGC);

            for iOri = 1:orientationsNum

                theConeModulationResponses = squeeze(theMRGCMosaicConeModulationsBasedResponses(iOri, :, :, theRGCindex));

                [theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 theConeModulationsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFspectra(...
                    theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                    theConeModulationResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    allowNonZeroMeanResponse, ...
                    false);

                % Allow zero mean of the pCurrent responses before fitting the sinusoid
                thePCurrentResponses = squeeze(theMRGCMosaicPhotocurrentBasedResponses(iOri, :, :, theRGCindex));
                theMeanPhotoCurrentResponsesAcrossEachTimeCourse = mean(thePCurrentResponses,2);
                thePCurrentResponses = bsxfun(@minus, thePCurrentResponses, theMeanPhotoCurrentResponsesAcrossEachTimeCourse);

                [thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFspectra(...
                    theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                    thePCurrentResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    allowNonZeroMeanResponse, ...
                    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses);

                coneModulationsSTF = squeeze(theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));
                photocurrentsSTF = squeeze(thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));

                theConeModulationsBasedBPIs(iRGC, iOri) = ...
                    bandPassIndex(stimParams.spatialFrequencyCPD, coneModulationsSTF);

                thePhotocurrentsBasedBPIs(iRGC, iOri) = ...
                    bandPassIndex(stimParams.spatialFrequencyCPD, photocurrentsSTF);

            end % iOri
        end % for iRGC
    end


    save(analyzedSTFsFileName, ...
        'stimParams', ...
        'targetRGCindices', ...
        'theCenterConePurities', ...
        'theCenterConeNumerosities', ...
        'theCenterConeDominances', ...
        'theSurroundConePurities', ...
        'theConeModulationsBasedTargetSTFamplitudeSpectra', ...
        'thePhotocurrentsBasedTargetSTFamplitudeSpectra', ...
        'theConeModulationsBasedTargetSTFphaseSpectra', ...
        'thePhotocurrentsBasedTargetSTFphaseSpectra', ...
        'theConeModulationsBasedBPIs', ...
        'thePhotocurrentsBasedBPIs');

    fprintf('\n\nSaved STFs to %s\n\n', analyzedSTFsFileName);



    if (visualizeConeExcitationVsPhotocurrentSTFs)

        for iRGC = 1:numel(targetRGCindices)

            theRGCindex = targetRGCindices(iRGC);

            if (all(stimParams.coneContrasts == [1 0 0]) && ...
               (theCenterConeDominances(iRGC) ~= cMosaic.LCONE_ID))

                % Skip. L-cone isolating stimulus but not L-cone dominated RF center
                continue;
            end

            if (all(stimParams.coneContrasts == [0 1 0]) && ...
               (theCenterConeDominances(iRGC) ~= cMosaic.MCONE_ID))

                % Skip. M-cone isolating stimulus but not M-cone dominated RF center
                continue;
            end

    
            % Visualize cone modulations based and photocurrents based STFs
            % and responses
            RGCMosaicAnalyzer.visualize.coneExcitationsVsPhotocurrentsSTF(...
                squeeze(theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC,:,:)), ...
                squeeze(thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC,:,:)), ...
                squeeze(theConeModulationsBasedTargetSTFphaseSpectra(iRGC,:,:)), ...
                squeeze(thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC,:,:)), ...
                squeeze(theMRGCMosaicConeModulationsBasedResponses(:, :, :, theRGCindex)), ...
                squeeze(theMRGCMosaicPhotocurrentBasedResponses(:, :, :, theRGCindex)), ...
                theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                stimParams, theRGCindex, ...
                max(abs(theMRGCMosaicConeModulationsBasedResponses(:))), ...
                max(abs(theMRGCMosaicPhotocurrentBasedResponses(:))), ...
                'exportPDFdirectory', exportPDFdirectory);


            if (1==2)
            
            hFig = figure(2); clf;
            set(hFig, 'Position', [1150 10 1100 1100], 'Color', [1 1 1]);


            for iOri = 1:orientationsNum
                coneModulationsSTF = squeeze(theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));
                photocurrentsSTF = squeeze(thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :));

                ax = subplot(1,orientationsNum+1,iOri);
                [p1,peakSFindex1] = max(coneModulationsSTF);
                [p2,peakSFindex2] = max(photocurrentsSTF);
                if (p1 > p2)
                    peakSFindex = peakSFindex1;
                else
                    peakSFindex = peakSFindex2;
                end

                factorToMatchLowSFs = coneModulationsSTF(peakSFindex)/photocurrentsSTF(peakSFindex);
                plot(ax, stimParams.spatialFrequencyCPD, coneModulationsSTF, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12);
                hold(ax, 'on');
                plot(ax, stimParams.spatialFrequencyCPD, photocurrentsSTF*factorToMatchLowSFs, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 12);
                set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], 'FontSize', 20, 'YLim', [0 1.0], 'YTick', 0:0.1:2);
                legend(ax, {'cone coneModulations based', 'photocurrent based'});

                if (theCenterConeDominances(iRGC) == cMosaic.LCONE_ID)
                    title(ax,sprintf('%d-(L-dominated) cone RF center\nstim. orientation = %d (degs)\n BPI(cone exc.):%2.2f, BPI(pCurrents):%2.2f', ...
                        theCenterConeNumerosities(iRGC), stimParams.orientationDegs(iOri), ...
                        theConeModulationsBasedBPIs(iRGC, iOri), thePhotocurrentsBasedBPIs(iRGC, iOri)));

                elseif (theCenterConeDominances(iRGC) == cMosaic.MCONE_ID)
                    title(ax,sprintf('%d-(M-dominated) cone RF center\nstim. orientation = %d (degs)\n BPI(cone exc.):%2.2f, BPI(pCurrents):%2.2f', ...
                        theCenterConeNumerosities(iRGC), stimParams.orientationDegs(iOri), ...
                        theConeModulationsBasedBPIs(iRGC, iOri), thePhotocurrentsBasedBPIs(iRGC, iOri)));

                elseif (theCenterConeDominances(iRGC) == cMosaic.SCONE_ID)
                    title(ax,sprintf('%d-(S-dominated) cone RF center\nstim. orientation = %d (degs)\n BPI(cone exc.):%2.2f, BPI(pCurrents):%2.2f', ...
                        theCenterConeNumerosities(iRGC), stimParams.orientationDegs(iOri), ...
                        theConeModulationsBasedBPIs(iRGC, iOri), thePhotocurrentsBasedBPIs(iRGC, iOri)));
                end

                grid(ax, 'on')
            end % for iOri

            drawnow;
            disp('Hit enter to continue')
            pause
            end


        end % for iRGC
    end % if (visualizeConeExcitationVsPhotocurrentSTFs)
end





function [theSTFamplitudeSpectrum, theSTFphaseDegsSpectrum, allFittedSFresponseTimeSeries] = computeSTFspectra(...
    temporalSupportSeconds, allSFresponseTimeSeries, temporalFrequencyHz, ...
    allowNonZeroMeanResponse, visualizeSinusoidalFits)


    theSFsNum = size(allSFresponseTimeSeries,1);
    theSTFamplitudeSpectrum = zeros(1, theSFsNum);
    theSTFphaseDegsSpectrum = zeros(1, theSFsNum);

    if (visualizeSinusoidalFits)
        hFig = figure(100); clf;
        set(hFig, 'Position', [10 10 1800 1150])
    end

    allFittedSFresponseTimeSeries = 0*allSFresponseTimeSeries;
    for iSF = 1:theSFsNum
        theResponseTimeSeries = allSFresponseTimeSeries(iSF,:);
        [theFittedSinusoid, fittedParams] = ...
            RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            temporalSupportSeconds, theResponseTimeSeries, temporalFrequencyHz, ...
            temporalSupportSeconds, ...
            'allowOffset', allowNonZeroMeanResponse);

        allFittedSFresponseTimeSeries(iSF,:) = theFittedSinusoid;
        theSTFamplitudeSpectrum(iSF) = fittedParams(1);
        theSTFphaseDegsSpectrum(iSF) = fittedParams(2);

        if (visualizeSinusoidalFits)
            ax = subplot(3, 7, iSF);
            plot(ax,temporalSupportSeconds, theResponseTimeSeries, 'ro', 'MarkerFaceColor', [1 0.5 0.5])
            hold(ax, 'on');
            plot(ax, temporalSupportSeconds, theFittedSinusoid, 'k-', 'LineWidth', 1.0);
            title(ax, sprintf('reponse amplitude: %2.2f\nsinusoid amplitude: %2.2f', max(abs(theResponseTimeSeries)), theSTFamplitudeSpectrum(iSF)));
            legend(ax, {'response', 'sinusoidal fit'});
            drawnow
        end

    end % iSF

    if (visualizeSinusoidalFits)
        disp('Hit enter to continue...')
        pause;
    end

end


function theBPI = bandPassIndex(theSpatialFrequencySupport, theSTF)

    [~, lowestSFindex] = min(theSpatialFrequencySupport(:));
    Ro = theSTF(lowestSFindex);
    [Rmax, indexOfMaxResponse] = max(theSTF(:));

    theBPI = Ro/Rmax;

    if (1==2)
        figure(9);
        plot(theSpatialFrequencySupport, theSTF, 'k-', 'LineWidth', 1.5);
        hold on;
        plot(theSpatialFrequencySupport(lowestSFindex), Ro, 'ro', 'MarkerSize', 16);
        plot(theSpatialFrequencySupport(indexOfMaxResponse), Rmax,  'bo', 'MarkerSize', 16);
        set(gca, 'XScale', 'log', 'XLim', [0.01 150]);
        title(sprintf('BPI = %2.3f', theBPI));
        drawnow;
        pause
    end
end

