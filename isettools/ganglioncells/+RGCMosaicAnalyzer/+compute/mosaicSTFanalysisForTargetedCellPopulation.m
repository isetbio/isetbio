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
    visualizeConeExcitationVsPhotocurrentSTFs)

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
        targetedCenterPurityRange, ...
        targetedCenterConeNumerosityRange, ...
        targetedSurroundPurityRange, ...
        targetedRadialEccentricityRange, ...
        theAnalyzedSTFsFileName, ...
        visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses, ...
        visualizeConeExcitationVsPhotocurrentSTFs);


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
    targetedCenterPurityRange, ...
    targetedCenterConeNumerosityRange, ...
    targetedSurroundPurityRange, ...
    targetedRadialEccentricityRange, ...
    analyzedSTFsFileName, ...
    visualizeSinusoidalFitsForPhotocurrentBasedMRGCresponses, ...
    visualizeConeExcitationVsPhotocurrentSTFs)

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
                 theConeModulationsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFamplitude(...
                    theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                    theConeModulationResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    false);

                [min(theConeModulationResponses(:)) max(theConeModulationResponses(:))]
                 % Allow zero mean of the pCurrent responses before fitting the sinusoid
                thePCurrentResponses = squeeze(theMRGCMosaicPhotocurrentBasedResponses(iOri, :, :, theRGCindex));
                [min(thePCurrentResponses(:)) max(thePCurrentResponses(:))]

                theMeanPhotoCurrentResponsesAcrossEachTimeCourse = mean(thePCurrentResponses,2)
                pause
                thePCurrentResponses = bsxfun(@minus, thePCurrentResponses, theMeanPhotoCurrentResponsesAcrossEachTimeCourse);

                [thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFamplitude(...
                    theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                    thePCurrentResponses, ...
                    stimParams.temporalFrequencyHz, ...
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
                 theConeModulationsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFamplitude(...
                    theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                    theConeModulationResponses, ...
                    stimParams.temporalFrequencyHz, ...
                    false);

                % Allow zero mean of the pCurrent responses before fitting the sinusoid
                thePCurrentResponses = squeeze(theMRGCMosaicPhotocurrentBasedResponses(iOri, :, :, theRGCindex));
                theMeanPhotoCurrentResponsesAcrossEachTimeCourse = mean(thePCurrentResponses,2);
                thePCurrentResponses = bsxfun(@minus, thePCurrentResponses, theMeanPhotoCurrentResponsesAcrossEachTimeCourse);

                [thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iOri, :), ...
                 thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC, iOri, :)] = computeSTFamplitude(...
                    theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                    thePCurrentResponses, ...
                    stimParams.temporalFrequencyHz, ...
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


            
            x = [];
            y = [];
            coneModulationsSTFMatrix = [];
            photocurrentsSTFMatrix = [];

            for iSF = 1:numel(stimParams.spatialFrequencyCPD)
                for iORI = 1:numel(stimParams.orientationDegs)
                    x(numel(x)+1) = iSF * cosd(stimParams.orientationDegs(iORI));
                    y(numel(y)+1) = iSF * sind(stimParams.orientationDegs(iORI));
                    coneModulationsSTFMatrix(numel(coneModulationsSTFMatrix)+1) = theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC,iORI,iSF);
                    photocurrentsSTFMatrix(numel(photocurrentsSTFMatrix)+1) = thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iORI, iSF);
                    x(numel(x)+1) = -x(numel(x));
                    y(numel(y)+1) = -y(numel(y));
                    coneModulationsSTFMatrix(numel(coneModulationsSTFMatrix)+1) = theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC,iORI,iSF);
                    photocurrentsSTFMatrix(numel(photocurrentsSTFMatrix)+1) = thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC, iORI, iSF);
                end
            end

            FconeModulations = scatteredInterpolant(x(:), y(:), coneModulationsSTFMatrix(:),'natural');
            Fphotocurrents = scatteredInterpolant(x(:), y(:), photocurrentsSTFMatrix(:),'natural');

            xx = -numel(stimParams.spatialFrequencyCPD):1:numel(stimParams.spatialFrequencyCPD);
            sfTicks = zeros(1, numel(xx));
            sfTickLabels = cell(1, numel(xx));
            for ix = 1:numel(xx)
                sfTicks(ix) = xx(ix);
                if (xx(ix) == 0)
                    sfTickLabels{ix} = sprintf('0');
                elseif (xx(ix)>=1)
                    sfTickLabels{ix} = sprintf('%2.0f', stimParams.spatialFrequencyCPD(xx(ix)));
                else
                    sfTickLabels{ix} = sprintf('%2.0f', -stimParams.spatialFrequencyCPD(abs(xx(ix))));
                end

            end

            yy = xx;
            [X,Y] = meshgrid(xx,yy);

            coneModulationsSTF2D = FconeModulations(X(:),Y(:));
            coneModulationsSTF2D = reshape(coneModulationsSTF2D, numel(xx), numel(yy));

            photocurrentsSTF2D = Fphotocurrents(X(:),Y(:));
            photocurrentsSTF2D = reshape(photocurrentsSTF2D, numel(xx), numel(yy));

            STFmax = max([max(coneModulationsSTFMatrix(:)) max(photocurrentsSTFMatrix(:))]);

            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 1750 800], 'Color', [1 1 1]);


    
            ax = subplot(2,4,3);
            hold(ax,'on')
            theLegends = cell(1,numel(stimParams.orientationDegs));
            pHandles = [];
            oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');
            for iORI = 1:numel(stimParams.orientationDegs)
                plot(stimParams.spatialFrequencyCPD, squeeze(theConeModulationsBasedTargetSTFamplitudeSpectra(iRGC,iORI,:)), ...
                    '-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5);
                theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
            end
            legend(ax, theLegends, 'Location', 'NorthEast');
            set(ax, 'FontSize', 16)
            xlabel('spatial frequency (c/deg)');
            ylabel('STF');
            title(ax, 'cone modulations');

            ax = subplot(2,4,4);
            hold(ax,'on')
            for iORI = 1:numel(stimParams.orientationDegs)
                plot(stimParams.spatialFrequencyCPD, squeeze(thePhotocurrentsBasedTargetSTFamplitudeSpectra(iRGC,iORI,:)), ...
                        '-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5);
            end
            legend(ax, theLegends, 'Location', 'NorthEast');
            set(ax, 'FontSize', 16)
            xlabel('spatial frequency (c/deg)');
            ylabel('STF');
            title(ax, 'photocurrents');

            [theMRGCMosaicConeModulationsResponseTemporalSupportSeconds(1) theMRGCMosaicConeModulationsResponseTemporalSupportSeconds(end)]
            [theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds(1) theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds(end)]
            size(theMRGCMosaicConeModulationsBasedResponses)
            size(theMRGCMosaicPhotocurrentBasedResponses)
            pause
            

            % The cone modulation response time-series
            ax = subplot(2,4,[1 5]);
            theYLims = max(abs(theMRGCMosaicConeModulationsBasedResponses(:))) * [-1 1];
            hold (ax, 'on');
            for iORI = 1:numel(stimParams.orientationDegs)
                if (iORI == 1)
                    sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
                else
                    sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
                end
                for iSF = 1:numel(stimParams.spatialFrequencyCPD)
                    theConeModulationsBasedResponse = squeeze(theMRGCMosaicConeModulationsBasedResponses(iORI, iSF, :, iRGC));

                    theConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse,...
                        theConeModulationsBasedTargetSTFphaseSpectra(iRGC, iORI, iSF), ...
                        theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, ...
                        stimParams, 1);
                    plot(ax, theMRGCMosaicConeModulationsResponseTemporalSupportSeconds, theConeModulationsBasedResponse, ...
                        '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
                end
            end
            grid(ax, 'on')
            set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', [0 theMRGCMosaicConeModulationsResponseTemporalSupportSeconds(end)]);
            xlabel('time (seconds)');
            ylabel('mRGC response');
            title(ax, 'cone modulations');

            ax = subplot(2,4,[2 6]);
            theYLims = max(abs(theMRGCMosaicPhotocurrentBasedResponses(:))) * [-1 1];
            hold (ax, 'on');
            for iORI = 1:numel(stimParams.orientationDegs)
                if (iORI == 1)
                    sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
                else
                    sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
                end
                for iSF = 1:numel(stimParams.spatialFrequencyCPD)
                    thePhotocurrentBasedResponse = squeeze(theMRGCMosaicPhotocurrentBasedResponses(iORI, iSF, :, iRGC));
                    thePhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse,...
                        thePhotocurrentsBasedTargetSTFphaseSpectra(iRGC, iORI, iSF), ...
                        theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, ...
                        stimParams, 0);
                    plot(ax, theMRGCMosaicPhotocurrentResponseTemporalSupportSeconds, thePhotocurrentBasedResponse, ...
                        '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
                end
            end
            grid(ax, 'on')
            set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', [0 theMRGCMosaicConeModulationsResponseTemporalSupportSeconds(end)]);
            xlabel('time (seconds)');
            ylabel('mRGC response');
            title(ax, 'photocurrents');      

            ax = subplot(2,4,7);
            imagesc(ax, xx,yy,coneModulationsSTF2D);
            set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
                    'YTick', sfTicks, 'YTickLabel', sfTickLabels);
            axis(ax, 'image');
            set(ax, 'CLim', [0 STFmax]);
            set(ax, 'FontSize', 16)
            title(ax, 'coneModulations');
            xlabel('spatial frequency (c/deg)');
            ylabel('spatial frequency (c/deg)');

            ax = subplot(2,4,8);
            imagesc(ax, xx,yy,photocurrentsSTF2D);
            set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
                    'YTick', sfTicks', 'YTickLabel', sfTickLabels);
            axis(ax, 'image');
            set(ax, 'CLim', [0 STFmax]);
            set(ax, 'FontSize', 16)
            title(ax, 'photocurrents')
            xlabel('spatial frequency (c/deg)');
            ylabel('spatial frequency (c/deg)');


            theFilename = sprintf('RGC_%d_nominalC_%2.0f%%_%2.0fCDM2_%2.1fHz.pdf', ...
                iRGC, stimParams.contrast*100, stimParams.backgroundLuminanceCdM2, stimParams.temporalFrequencyHz);
            
            NicePlot.exportFigToPDF(theFilename, hFig, 300);



            hFig = figure(2); clf;
            set(hFig, 'Position', [10 10 1100 1100], 'Color', [1 1 1]);


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
                        theConeconeModulationsBasedBPIs(iRGC, iOri), thePhotocurrentsBasedBPIs(iRGC, iOri)));

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
        end % for iRGC
    end % if (visualizeConeExcitationVsPhotocurrentSTFs)
end


function theResponse = phaseAlignResponse(theResponse, theResponsePhaseDegs, theTemporalSupportSeconds, stimParams, extraSamples)

    sizeResponse = size(theResponse);
    thePeriodSeconds = 1/stimParams.temporalFrequencyHz;
    theStepPhaseDegs = 360*theTemporalSupportSeconds/thePeriodSeconds;
    [~,theShiftAmount] = min(abs(theStepPhaseDegs-theResponsePhaseDegs));
    if (extraSamples > 0)
        theResponse = theResponse(1:end-extraSamples);
    end
    theResponse  = circshift(theResponse , -theShiftAmount);

    if (extraSamples > 0)
        theResponse = [theResponse(:);  theResponse(1:extraSamples)];
    end

    theResponse = reshape(theResponse, sizeResponse);
end

function [theSTFamplitudeSpectrum, theSTFphaseDegsSpectrum] = computeSTFamplitude(...
    temporalSupportSeconds, allSFresponseTimeSeries, temporalFrequencyHz, ...
    visualizeSinusoidalFits)

    theSFsNum = size(allSFresponseTimeSeries,1);
    theSTFamplitude = zeros(1, theSFsNum);
    theSTFphaseDegs = zeros(1, theSFsNum);

    if (visualizeSinusoidalFits)
        hFig = figure(100); clf;
        set(hFig, 'Position', [10 10 1800 1150])
    end

    for iSF = 1:theSFsNum
        theResponseTimeSeries = allSFresponseTimeSeries(iSF,:);
        [theFittedSinusoid, fittedParams] = ...
            RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            temporalSupportSeconds, theResponseTimeSeries, temporalFrequencyHz, temporalSupportSeconds);

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

