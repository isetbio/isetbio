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

    temporalFrequenciesExamined = TTFparamsStruct.tfSupport;

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






%
%RGCMosaicAnalyzer.compute.MRGCtemporalFiltersFromPhotocurrentsBasedTTF
%
%
function MRGCtemporalFiltersFromPhotocurrentsBasedTTFnew(...
    stimulusShape, ...
    theTargetRGCwithIndex, ...
    theMRGCMosaicTTFResponsesFullFileName, ...
    theAnalyzedTTFsFullFileName, ...
    visualizeMosaicResponses)


    load(theMRGCMosaicTTFResponsesFullFileName, ...
        'theMRGCMosaic', 'stimParams', 'TTFparamsStruct', ...
        'computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex', ...
        'theMRGCMosaicTTFresponsesAllConditions', ...
        'theMRGCMosaicTemporalSupportSecondsAllConditions');

    temporalFrequenciesExamined = TTFparamsStruct.tfSupport;
    temporalFrequenciesExamined = temporalFrequenciesExamined(temporalFrequenciesExamined< 120);

    responseRange = max(max(abs(squeeze(theMRGCMosaicTTFresponsesAllConditions(:,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex)))))*[-1 1];

    theTTFamplitude = zeros(1, numel(temporalFrequenciesExamined));
    theTTFphaseDegs = zeros(1, numel(temporalFrequenciesExamined));

    for iTF = 1:numel(temporalFrequenciesExamined)

        % Retrieve the photocurrent response data for this TF
        theTemporalSupportSecondsForThisTF = squeeze(theMRGCMosaicTemporalSupportSecondsAllConditions(iTF,:));
        theMRGCMosaicPhotocurrentResponsesForThisTF = squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex));

        nanIndices = find(isnan(squeeze(theMRGCMosaicTTFresponsesAllConditions(iTF,:,computePhotocurrentResponsesOnlyForInputsToSingleRGCwithIndex))));
        if (isempty(nanIndices))
            theTimeBins = 1:size(theMRGCMosaicTTFresponsesAllConditions,2);
        else
            theTimeBins = 1:(nanIndices(1)-1);
        end

        theTemporalSupportSecondsForThisTF = theTemporalSupportSecondsForThisTF(theTimeBins);
        theMRGCMosaicPhotocurrentResponsesForThisTF = theMRGCMosaicPhotocurrentResponsesForThisTF(theTimeBins);

        assert(isempty(find(isnan(theMRGCMosaicPhotocurrentResponsesForThisTF(:)))), 'did not remove all nan part of the response');
        assert(isempty(find(isnan(theTemporalSupportSecondsForThisTF(:)))), 'did not remove all nan part of the response temporal support');

        dt = theTemporalSupportSecondsForThisTF(2)-theTemporalSupportSecondsForThisTF(1);
        temporalSupport1Msec = 0:1/1000:(theTemporalSupportSecondsForThisTF(end)+dt);

        [theFittedResponse, fittedParams] = RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries(...
            theTemporalSupportSecondsForThisTF, ...
            theMRGCMosaicPhotocurrentResponsesForThisTF, ...
            temporalFrequenciesExamined(iTF), ...
            temporalSupport1Msec, ...
            'allowOffset', true);

        theTTFamplitude(iTF) = fittedParams(1);
        theTTFphaseDegs(iTF) = fittedParams(2);

        if (visualizeMosaicResponses)
            hFig = figure(1); clf;
            ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
            theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
            ax = theAxes{1,1};

            if (theTemporalSupportSecondsForThisTF(end) > 1)
                dTick = 0.4;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.5)
                dTick = 0.2;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.3)
                dTick = 0.1;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.1)
                dTick = 0.05;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.06)
                dTick = 0.02;
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.03)
                dTick = 0.01;
            else
                dTick = 0.005;
            end
            xTicks = 0:dTick:2.0;
            if (theTemporalSupportSecondsForThisTF(end) > 0.5)
                xTickLabels = sprintf('%.1f\n', xTicks);
            elseif (theTemporalSupportSecondsForThisTF(end) > 0.1)
                xTickLabels = sprintf('%.2f\n', xTicks);
            else
                xTickLabels = sprintf('%.3f\n', xTicks);
            end

            XLims = [0 theTemporalSupportSecondsForThisTF(end)];
            YLims = [-2 2];

            plot(ax,theTemporalSupportSecondsForThisTF, theMRGCMosaicPhotocurrentResponsesForThisTF, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
            hold(ax, 'on')
            plot(ax, temporalSupport1Msec, theFittedResponse, 'k-', 'LineWidth', 1.5);
            set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', xTicks, 'YTick', [-20:1:2], 'XTickLabel', xTickLabels);
            grid(ax, 'on');
            xlabel(ax,'time (sec)');
            ylabel(ax,'pAmps');
            title(ax, sprintf('%s, TF = %2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF)));

            % Finalize figure using the Publication-Ready format
            PublicationReadyPlotLib.applyFormat(ax,ff);
            PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

            % Export
            visualizationPDFfileName = sprintf('%s_TF_%2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF));
            exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
            pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
            theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

            % Generate the path if we need to
            RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                pdfExportRootDir, theVisualizationPDFfilename, ...
                'generateMissingSubDirs', true);

            thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
            NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');
        end % visualizeMosaicResponses


        if (1==2)
            hFig = figure(2); clf
            set(hFig, 'Position', [10 10 1500 500], 'Name', sprintf('Responses to %s stimulus', stimulusShape));
            axResponse = subplot(1,3,1);
            axTTFamplitude = subplot(1,3,2);
            axTTFphase = subplot(1,3,3);

            plot(axResponse,theTemporalSupportSecondsForThisTF*1e3, theMRGCMosaicPhotocurrentResponsesForThisTF, 'k.');
            hold(axResponse, 'on')
            plot(axResponse, temporalSupport1Msec*1e3, theFittedResponse, 'r-');
            set(axResponse, 'XLim', [0 theTemporalSupportSecondsForThisTF(end)]*1e3, 'YLim', max(abs(theMRGCMosaicPhotocurrentResponsesForThisTF(:))) * [-1 1]);
            xlabel(axResponse,'time (msec)');
            ylabel(axResponse,'pAmps');
            title(axResponse, sprintf('%s, TF = %2.2fHz', stimulusShape, temporalFrequenciesExamined(iTF)));

            plot(axTTFamplitude, temporalFrequenciesExamined, theTTFamplitude, 'ko-');
            set(axTTFamplitude, 'XLim', [0.3 300], 'XTick', [0.3 1 3 10 30 100 300], 'XScale', 'log');
            grid(axTTFamplitude, 'on');
            xlabel(axTTFamplitude,'frequency (Hz)');
            ylabel(axTTFamplitude,'amplitude');

            plot(axTTFphase, temporalFrequenciesExamined, theTTFphaseDegs, 'ko-');
            set(axTTFphase, 'XLim', [0.3 300], 'XTick', [0.3 1 3 10 30 100 300], 'XScale', 'log');
            grid(axTTFphase, 'on');
            xlabel(axTTFphase,'frequency (Hz)');
            ylabel(axTTFphase,'phase (degs)');

            drawnow;
        end

    end % iTF

    % Massage the phase
    theTTFphaseDegs = -theTTFphaseDegs;

    % Phase in radians
    theTTFphaseRadians = theTTFphaseDegs/180*pi;

    % Complex TTF from amplitude and phase
    thePhotocurrentBasedMRGCcellTTF = theTTFamplitude .* exp(1i * theTTFphaseRadians);



    % Generate a typical Benardete&Kaplan mRGC center TTFs
    nL_tL_product = 48;


    % Center TTF params (ON-center cell data from Figure 6)
    centerParamsFig6(1) = 184.2;        % gain (A)
    centerParamsFig6(2) = 0.69;         % high-pass gain (Hs)
    centerParamsFig6(3) = 18.61;        % high-pass time constant (msec) (Ts)
    centerParamsFig6(4) = 38;           % low-pass stages num (Nl)
    centerParamsFig6(5) = nL_tL_product/centerParamsFig6(4);         % low-pass time constant (msec) (Tl)
    centerParamsFig6(6) = 4.0;          % delay (msec) (D)
    theBenardeteKaplanFigure6MRGCcenterTTF = highPassNstageLowPassTTF(centerParamsFig6, temporalFrequenciesExamined);
    theBenardeteKaplanFigure6MRGCcenterTTFamplitude = abs(theBenardeteKaplanFigure6MRGCcenterTTF);

    % Values from Figure 7
    centerParamsFig7(1) = 184.2;        % gain (A)
    centerParamsFig7(2) = 0.70;         % high-pass gain (Hs)
    centerParamsFig7(3) = 37.30;        % high-pass time constant (msec) (Ts)
    centerParamsFig7(4) = 32;           % low-pass stages num (Nl)
    centerParamsFig7(5) = nL_tL_product/centerParamsFig7(4);         % low-pass time constant (msec) (Tl)
    centerParamsFig7(6) = 4.0;          % delay (msec) (D)
    theBenardeteKaplanFigure7MRGCcenterTTF = highPassNstageLowPassTTF(centerParamsFig7, temporalFrequenciesExamined);
    theBenardeteKaplanFigure7MRGCcenterTTFamplitude = abs(theBenardeteKaplanFigure7MRGCcenterTTF);

    % Surround TTF params (ON-center cell data from Figure 6)
    surroundParamsFig6(1) = -125.33;       % gain (A)
    surroundParamsFig6(2) = 0.56;         % high-pass gain (Hs)
    surroundParamsFig6(3) = 33.28;        % high-pass time constant (msec) (Ts)
    surroundParamsFig6(4) = 124;          % low-pass stages num (Nl)
    surroundParamsFig6(5) = nL_tL_product/surroundParamsFig6(4);         % low-pass time constant (msec) (Tl)
    surroundParamsFig6(6) = 4.03;         % delay (msec) (D)
    theBenardeteKaplanFigure6MRGCsurroundTTF = highPassNstageLowPassTTF(surroundParamsFig6, temporalFrequenciesExamined);
    theBenardeteKaplanFigure6MRGCsurroundTTFamplitude = abs(theBenardeteKaplanFigure6MRGCsurroundTTF);

    % Values from Figure 7
    surroundParamsFig7(1) = -centerParamsFig7(1)*0.64;       % gain (A)
    surroundParamsFig7(2) = 0.41;         % high-pass gain (Hs)
    surroundParamsFig7(3) = 42.49;        % high-pass time constant (msec) (Ts)
    surroundParamsFig7(4) = 93;          % low-pass stages num (Nl)
    surroundParamsFig7(5) = nL_tL_product/surroundParamsFig7(4);         % low-pass time constant (msec) (Tl)
    surroundParamsFig7(6) = 4.0;          % delay (msec) (D)
    theBenardeteKaplanFigure7MRGCsurroundTTF = highPassNstageLowPassTTF(surroundParamsFig7, temporalFrequenciesExamined);
    theBenardeteKaplanFigure7MRGCsurroundTTFamplitude = abs(theBenardeteKaplanFigure7MRGCsurroundTTF);



    % photocurrents-based TIR
    zeroPaddingLength = 512;
    performFFTshift = false;
    delayMilliseconds = 20;
    thePhotocurrentImpulseResponseData = temporalTransferFunctionToImpulseResponseFunction(...
        thePhotocurrentBasedMRGCcellTTF, temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);

    % Benardete & Kaplan Fig 6&7 TIR (center)
    delayMilliseconds = 0;
    theBenardeteKaplanFigure6MRGCcenterTIR  = temporalTransferFunctionToImpulseResponseFunction(...
        theBenardeteKaplanFigure6MRGCcenterTTF , temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);

    delayMilliseconds = 0;
    theBenardeteKaplanFigure7MRGCcenterTIR  = temporalTransferFunctionToImpulseResponseFunction(...
        theBenardeteKaplanFigure7MRGCcenterTTF , temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);


    % Benardete & Kaplan Fig 6&7 TIR (surrounds)
    delayMilliseconds = 0;
    theBenardeteKaplanFigure6MRGCsurroundTIR  = temporalTransferFunctionToImpulseResponseFunction(...
        theBenardeteKaplanFigure6MRGCsurroundTTF , temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);

    delayMilliseconds = 0;
    theBenardeteKaplanFigure7MRGCsurroundTIR  = temporalTransferFunctionToImpulseResponseFunction(...
        theBenardeteKaplanFigure7MRGCsurroundTTF , temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);


    % Derive inner retina TTF
    deconvolutionMode = 'amplitudeDivisionOnly';
    deconvolutionMode = 'fullDivision';

    mRGCtoMatch = 'fromFig6';
    %mRGCtoMatch = 'fromFig7';

    switch (deconvolutionMode)
        case 'amplitudeDivisionOnly'
            % Deconvolve by dividing with the amplitude
            if (strcmp(stimulusShape, 'spot'))
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        theInnerRetinaTTFamplitude = theBenardeteKaplanFigure6MRGCcenterTTFamplitude ./ abs(thePhotocurrentBasedMRGCcellTTF);
                        theInnerRetinaTTFphaseRadians = angle(theBenardeteKaplanFigure6MRGCcenterTTF);
                    case 'fromFig7'
                        theInnerRetinaTTFamplitude = theBenardeteKaplanFigure7MRGCcenterTTFamplitude ./ abs(thePhotocurrentBasedMRGCcellTTF);
                        theInnerRetinaTTFphaseRadians = angle(theBenardeteKaplanFigure7MRGCcenterTTF);
                    otherwise
                        error('Unkown mRGCtoMatch: ''%s''.', mRGCtoMatch)
                end % switch (mRGCtoMatch)

                % Complex TTF from amplitude and phase
                theInnerRetinaTTF = theInnerRetinaTTFamplitude .* exp(1i * theInnerRetinaTTFphaseRadians );
            else
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        theInnerRetinaTTFamplitude = theBenardeteKaplanFigure6MRGCsurroundTTFamplitude ./ abs(thePhotocurrentBasedMRGCcellTTF);
                        theInnerRetinaTTFphaseRadians = angle(theBenardeteKaplanFigure6MRGCsurroundTTF);
                    case 'fromFig7'
                        theInnerRetinaTTFamplitude = theBenardeteKaplanFigure7MRGCsurroundTTFamplitude ./ abs(thePhotocurrentBasedMRGCcellTTF);
                        theInnerRetinaTTFphaseRadians = angle(theBenardeteKaplanFigure7MRGCsurroundTTF);
                    otherwise
                        error('Unkown mRGCtoMatch: ''%s''.', mRGCtoMatch);
                end % switch (mRGCtoMatch)

                % Complex TTF from amplitude and phase
                theInnerRetinaTTF = theInnerRetinaTTFamplitude .* exp(1i * theInnerRetinaTTFphaseRadians );
            end

        case 'fullDivision'
            if (strcmp(stimulusShape, 'spot'))
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        theInnerRetinaTTF = theBenardeteKaplanFigure6MRGCcenterTTF ./ thePhotocurrentBasedMRGCcellTTF;
                        theInnerRetinaTTFamplitude = abs(theInnerRetinaTTF);
                    case 'fromFig7'
                        theInnerRetinaTTF = theBenardeteKaplanFigure7MRGCcenterTTF ./ thePhotocurrentBasedMRGCcellTTF;
                        theInnerRetinaTTFamplitude = abs(theInnerRetinaTTF);
                    otherwise
                        error('Unkown mRGCtoMatch: ''%s''.', mRGCtoMatch);
                end
            else
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        theInnerRetinaTTF = theBenardeteKaplanFigure6MRGCsurroundTTF ./ thePhotocurrentBasedMRGCcellTTF;
                        theInnerRetinaTTFamplitude = abs(theInnerRetinaTTF);
                    case 'fromFig7'
                        theInnerRetinaTTF = theBenardeteKaplanFigure7MRGCsurroundTTF ./ thePhotocurrentBasedMRGCcellTTF;
                        theInnerRetinaTTFamplitude = abs(theInnerRetinaTTF);
                    otherwise
                        error('Unkown mRGCtoMatch: ''%s''.', mRGCtoMatch);
                end
            end

        otherwise
            error('Unknown deconvolution mode: ''%s''.', deconvolutionMode)
    end


    % Derive the cascade TTF
    theCascadePhotocurrentsInnerRetinaTTF = thePhotocurrentBasedMRGCcellTTF .* theInnerRetinaTTF;
    theCascadePhotocurrentsInnerRetinaTTFamplitude = abs(theCascadePhotocurrentsInnerRetinaTTF);

    % Derive the inner retina TIR
    theInnerRetinaTIR = temporalTransferFunctionToImpulseResponseFunction(...
        theInnerRetinaTTF, temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);

    % Derive the cascale TIR
    delayMilliseconds = 0;
    theCascadePhotocurrentInnerRetinaTIR = temporalTransferFunctionToImpulseResponseFunction(...
        theCascadePhotocurrentsInnerRetinaTTF, temporalFrequenciesExamined, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift);


    % Build the figure in stages and export each stage separately
    for iBuildUpStage = 1:4

        allPlotHandles = [];
        allLegends = {};

        % Plot the TTF
        hFig = figure(60); clf;
        ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        XLims = [0.3 200];
        YLims = [0 1.01];

        p1 = plot(ax,temporalFrequenciesExamined, theTTFamplitude/max(theTTFamplitude), 'ro-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
        allPlotHandles(numel(allPlotHandles)+1) = p1;
        allLegends{numel(allLegends)+1} = 'photocurrents-derived';

        if (iBuildUpStage > 1)
            hold(ax, 'on');
            if (strcmp(stimulusShape, 'spot'))
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure6MRGCcenterTTFamplitude/max(theBenardeteKaplanFigure6MRGCcenterTTFamplitude), 'k-', ...
                            'LineWidth', 4, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        p2 = plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure6MRGCcenterTTFamplitude/max(theBenardeteKaplanFigure6MRGCcenterTTFamplitude), 'c-', ...
                            'LineWidth', 2, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        allPlotHandles(numel(allPlotHandles)+1) = p2;
                        allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 6)';
                    case 'fromFig7'
                        plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure7MRGCcenterTTFamplitude/max(theBenardeteKaplanFigure7MRGCcenterTTFamplitude), 'k-', ...
                            'LineWidth', 4, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        p2 = plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure7MRGCcenterTTFamplitude/max(theBenardeteKaplanFigure7MRGCcenterTTFamplitude), 'c-', ...
                            'LineWidth', 2, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        allPlotHandles(numel(allPlotHandles)+1) = p2;
                        allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 7)';
                    end
            else
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure6MRGCsurroundTTFamplitude/max(theBenardeteKaplanFigure6MRGCsurroundTTFamplitude), 'k-', ...
                            'LineWidth', 4, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        p2 = plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure6MRGCsurroundTTFamplitude/max(theBenardeteKaplanFigure6MRGCsurroundTTFamplitude), 'c-', ...
                            'LineWidth', 2, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        allPlotHandles(numel(allPlotHandles)+1) = p2;
                        allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 6)';
                    case 'fromFig7'
                        plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure7MRGCsurroundTTFamplitude/max(theBenardeteKaplanFigure7MRGCsurroundTTFamplitude), 'k-', ...
                            'LineWidth', 4, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        p2 = plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure7MRGCsurroundTTFamplitude/max(theBenardeteKaplanFigure7MRGCsurroundTTFamplitude), 'c-', ...
                            'LineWidth', 2, ...
                            'MarkerSize', 14, 'MarkerFaceColor', [1 0.75 0.75]);
                        allPlotHandles(numel(allPlotHandles)+1) = p2;
                        allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 7)';
                    end
            end



            if (iBuildUpStage > 2)
                %p3 = plot(ax,temporalFrequenciesExamined, theBenardeteKaplanFigure7MRGCcenterTTFamplitude/max(theBenardeteKaplanFigure7MRGCcenterTTFamplitude), 'k--', ...
                %    'LineWidth', 1.5, ...
                %    'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
                plot(ax,temporalFrequenciesExamined, theInnerRetinaTTFamplitude/max(theInnerRetinaTTFamplitude), 'k-', ...
                    'LineWidth', 4, ...
                    'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
                p4 = plot(ax,temporalFrequenciesExamined, theInnerRetinaTTFamplitude/max(theInnerRetinaTTFamplitude), 'y-', ...
                    'LineWidth', 2, ...
                    'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);

                allPlotHandles(numel(allPlotHandles)+1) = p4;
                allLegends{numel(allLegends)+1} = 'inner retina filter (derived)';

                if (iBuildUpStage > 3)
                    p5 =  plot(ax,temporalFrequenciesExamined, theCascadePhotocurrentsInnerRetinaTTFamplitude/max(theCascadePhotocurrentsInnerRetinaTTFamplitude), 'd', ...
                        'LineWidth', 2.0, ...
                        'MarkerSize', 10, 'MarkerEdgeColor', 0.3*[1 0.8 0.5], ...
                        'MarkerFaceColor', 'none');

                    allPlotHandles(numel(allPlotHandles)+1) = p5;
                    allLegends{numel(allLegends)+1} = 'photocurrent x inner retina cascade';

                end % if (iBuildUpStage > 3)
            end  % if (iBuildUpStage > 2)
        end % if (iBuildUpStage > 1)

        legend(ax, allPlotHandles, allLegends, 'Location', 'SouthWest', 'NumColumns', 1);
        xlabel(ax, 'frequency (Hz)')
        set(ax, 'XLim', XLims, 'XTick', [0.3 1 3 10 30 100], 'YLim', YLims, 'XScale', 'log');
        grid(ax, 'on')
        ff.legendBox = 'on';
        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
        %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

        % Export
        visualizationPDFfileName = sprintf('%s_photocurrentsBasedTTF_BuildUpStage_%d', stimulusShape, iBuildUpStage);
        exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
        pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

        thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
        NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');
    end  % for iBuildUpStage = 1:4



    addLegends = true;
    %addLegends = false;

    % Build the figure in stages and export each stage separately
    for iBuildUpStage = 1:4

        if (addLegends)
            allPlotHandles = [];
            allLegends = {};
        end

        % Plot TTF-derived impulse response function
        hFig = figure(61); clf;
        ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};

        XLims = [0 200]
        YLims = 1.05*[-1 1];

        p1 = plot(ax,thePhotocurrentImpulseResponseData.temporalSupportSeconds*1e3, ...
            thePhotocurrentImpulseResponseData.weights/max(abs(thePhotocurrentImpulseResponseData.weights)), 'ro-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);

        if (addLegends)
            allPlotHandles(numel(allPlotHandles)+1) = p1;
            allLegends{numel(allLegends)+1} = 'photocurrents-derived';
        end

        if (iBuildUpStage > 1)
            hold(ax, 'on');

            if (strcmp(stimulusShape, 'spot'))
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        plot(ax, theBenardeteKaplanFigure6MRGCcenterTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure6MRGCcenterTIR.weights/max(abs(theBenardeteKaplanFigure6MRGCcenterTIR.weights)), 'k-', ...
                            'LineWidth', 4);
                        p2 = plot(ax, theBenardeteKaplanFigure6MRGCcenterTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure6MRGCcenterTIR.weights/max(abs(theBenardeteKaplanFigure6MRGCcenterTIR.weights)), 'c-', ...
                            'LineWidth', 2.0);
                        if (addLegends)
                            allPlotHandles(numel(allPlotHandles)+1) = p2;
                            allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 6)';
                        end
                    case 'fromFig7'
                        plot(ax, theBenardeteKaplanFigure7MRGCcenterTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure7MRGCcenterTIR.weights/max(abs(theBenardeteKaplanFigure7MRGCcenterTIR.weights)), 'k-', ...
                            'LineWidth', 4);
                        p2 = plot(ax, theBenardeteKaplanFigure7MRGCcenterTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure7MRGCcenterTIR.weights/max(abs(theBenardeteKaplanFigure7MRGCcenterTIR.weights)), 'c-', ...
                            'LineWidth', 2.0);
                        if (addLegends)
                            allPlotHandles(numel(allPlotHandles)+1) = p2;
                            allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 7)';
                        end
                    end

            else
                switch (mRGCtoMatch)
                    case 'fromFig6'
                        plot(ax, theBenardeteKaplanFigure6MRGCsurroundTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure6MRGCsurroundTIR.weights/max(abs(theBenardeteKaplanFigure6MRGCsurroundTIR.weights)), 'k-', ...
                            'LineWidth', 4);
                        p2 = plot(ax, theBenardeteKaplanFigure6MRGCsurroundTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure6MRGCsurroundTIR.weights/max(abs(theBenardeteKaplanFigure6MRGCsurroundTIR.weights)), 'c-', ...
                            'LineWidth', 2.0);
                        if (addLegends)
                            allPlotHandles(numel(allPlotHandles)+1) = p2;
                            allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 6)';
                        end
                    case 'fromFig7'
                        plot(ax, theBenardeteKaplanFigure7MRGCsurroundTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure7MRGCsurroundTIR.weights/max(abs(theBenardeteKaplanFigure6MRGCsurroundTIR.weights)), 'k-', ...
                            'LineWidth', 4);
                        p2 = plot(ax, theBenardeteKaplanFigure7MRGCsurroundTIR.temporalSupportSeconds*1e3, ...
                            theBenardeteKaplanFigure7MRGCsurroundTIR.weights/max(abs(theBenardeteKaplanFigure7MRGCsurroundTIR.weights)), 'c-', ...
                            'LineWidth', 2.0);
                        if (addLegends)
                            allPlotHandles(numel(allPlotHandles)+1) = p2;
                            allLegends{numel(allLegends)+1} = 'Benardete&Kaplan (Fig 7)';
                        end
                    end
            end


            if (iBuildUpStage > 2)
                %p3 = plot(ax, theBenardeteKaplanFigure7MRGCcenterTIR.temporalSupportSeconds*1e3, ...
                %    theBenardeteKaplanFigure7MRGCcenterTIR.weights/max(abs(theBenardeteKaplanFigure7MRGCcenterTIR.weights)), 'k--', ...
                %    'LineWidth', 1.5);
                plot(ax, theInnerRetinaTIR.temporalSupportSeconds*1e3, ...
                    theInnerRetinaTIR.weights/max(abs(theInnerRetinaTIR.weights)), 'k-', ...
                    'LineWidth', 5);
                hold(ax, 'on')
                p4 = plot(ax, theInnerRetinaTIR.temporalSupportSeconds*1e3, ...
                    theInnerRetinaTIR.weights/max(abs(theInnerRetinaTIR.weights)), 'y-', ...
                    'LineWidth', 3);
                if (addLegends)
                    allPlotHandles(numel(allPlotHandles)+1) = p4;
                    allLegends{numel(allLegends)+1} = 'inner retina filter (derived)';
                end

                if (iBuildUpStage > 3)
                    p5 = plot(ax, theCascadePhotocurrentInnerRetinaTIR.temporalSupportSeconds*1e3, ...
                        theCascadePhotocurrentInnerRetinaTIR.weights/max(abs(theCascadePhotocurrentInnerRetinaTIR.weights)), 'd', ...
                        'LineWidth', 2.0, ...
                        'MarkerSize', 10, 'MarkerEdgeColor', 0.3*[1 0.8 0.5], ...
                        'MarkerFaceColor', 'none');

                    if (addLegends)
                        allPlotHandles(numel(allPlotHandles)+1) = p5;
                        allLegends{numel(allLegends)+1} = 'photocurrent x inner retina cascade';
                    end

                end % if (iBuildUpStage > 3)
            end  % if (iBuildUpStage > 2)
        end % if (iBuildUpStage > 1)

        if (addLegends)
            if (strcmp(stimulusShape, 'spot'))
                legend(ax, allPlotHandles, allLegends, 'Location', 'SouthEast', 'NumColumns', 1);
            else
                legend(ax, allPlotHandles, allLegends, 'Location', 'NorthEast', 'NumColumns', 1);
            end
        end

        xlabel(ax, 'time (msec)')
        set(ax, 'XLim', XLims, 'XTick', 0:50:200, 'YLim', YLims);
        grid(ax, 'on')
        ff.legendBox = 'on';

        % Finalize figure using the Publication-Ready format
        PublicationReadyPlotLib.applyFormat(ax,ff);
        %PublicationReadyPlotLib.offsetAxes(ax, ff, XLims, YLims);

        % Export
        visualizationPDFfileName = sprintf('%s_photocurrentsBasedIR_BuildUpStage_%d', stimulusShape, iBuildUpStage);
        exportVisualizationPDFdirectory = 'temporalResponseGenerationPDFs';
        pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
        theVisualizationPDFfilename = fullfile(exportVisualizationPDFdirectory, sprintf('%s.pdf', visualizationPDFfileName));

        thePDFfileName = fullfile(pdfExportRootDir, theVisualizationPDFfilename);
        NicePlot.exportFigToPDF(thePDFfileName, hFig, 300, 'beVerbose');

    end % for iBuildUpStage = 1:4
end



function previous()
    pause


    hold on;
    p2 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, theMacaqueCenterImpulseResponseData.weights/max(theMacaqueCenterImpulseResponseData.weights), 'r-');
    p3 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, theInnerRetinaImpulseResponseData.weights/max(theInnerRetinaImpulseResponseData.weights), 'b-');
    p4 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, prediction/max(prediction), 'ms-');
    legend([p1, p2, p3, p4], {'pCurrent', 'macaque', 'inner retina', 'prediction macaque'})
    xlabel('time (msec)')
    set(gca, 'XLim', [0 500], 'XTick', 0:10:1000);
    grid on




    nL_tL_product = 48;
    centerParams(1) = 184.2;        % gain (A)

    centerParams(2) = 0.69;         % high-pass gain (Hs)
    centerParams(3) = 18.61;        % high-pass time constant (msec) (Ts)

    centerParams(4) = 38;           % low-pass stages num (Nl)
    centerParams(5) = nL_tL_product/centerParams(4);         % low-pass time constant (msec) (Tl)

    centerParams(6) = 4.0;          % delay (msec) (D)

    theMacaqueCenterTTF = highPassNstageLowPassTTF(centerParams, temporalFrequenciesExamined);

    theMacaqueCenterImpulseResponseData = temporalTransferFunctionToImpulseResponseFunction(...
            theMacaqueCenterTTF, temporalFrequenciesExamined, false);

    [theInnerRetinaImpulseResponseData.weights,tmpMacaqueWeights] = ...
        deconv(theMacaqueCenterImpulseResponseData.weights, ...
            thePhotocurrentImpulseResponseData.weights, 'same', ...
            Method='least-squares', RegularizationFactor=0.005*4);

    % We see that witout regularization (Factor = 0), we get noise.
    % With regularization, we get better results, but prediction is not good
    % Predict the macaque IR
    prediction = conv(theInnerRetinaImpulseResponseData.weights, thePhotocurrentImpulseResponseData.weights, 'same');

    figure(22); clf
    plot(temporalFrequenciesExamined, unwrap(angle(theMacaqueCenterTTF))/pi*180, 'ko-')
    hold on;
    plot(temporalFrequenciesExamined, unwrap(angle(thePhotocurrentBasedMRGCcellTTF))/pi*180, 'r.-')
    plot(temporalFrequenciesExamined, unwrap(theTTFphaseRadians)/pi*180, 'r--')





    % Plot impulse response functions
    figure(60); clf;
    p1 = plot(thePhotocurrentImpulseResponseData.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseData.weights/max(thePhotocurrentImpulseResponseData.weights), 'ko-');
    hold on;
    p2 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, theMacaqueCenterImpulseResponseData.weights/max(theMacaqueCenterImpulseResponseData.weights), 'r-');
    p3 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, theInnerRetinaImpulseResponseData.weights/max(theInnerRetinaImpulseResponseData.weights), 'b-');
    p4 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, prediction/max(prediction), 'ms-');
    legend([p1, p2, p3, p4], {'pCurrent', 'macaque', 'inner retina', 'prediction macaque'})
    xlabel('time (msec)')
    set(gca, 'XLim', [0 500], 'XTick', 0:10:1000);
    grid on



    % Lets try a high pass filter
    highPassParams(1) = 1.0 % gain (A)
    highPassParams(2) = 100;            % high-pass stages num
    highPassParams(3) = 100;          % high-pass time constant (msec) (Ts)
    highPassParams(4) = 50.0;          % delay (msec) (D)

    % Lets try a high pass filter TTF
    theHighPassTTF = highPassNstageTTF(highPassParams, temporalFrequenciesExamined)

    theHighPassImpulseResponseData = ...
        temporalTransferFunctionToImpulseResponseFunction(theHighPassTTF, temporalFrequenciesExamined, false);

    % We see that witout regularization (Factor = 0), we get noise.
    % With regularization, we get better results, but prediction is not good
    % Predict the macaque IR
    highPassPrediction = conv(theHighPassImpulseResponseData.weights, thePhotocurrentImpulseResponseData.weights, 'same');



     % Plot impulse response functions
    figure(61); clf;
    p1 = plot(thePhotocurrentImpulseResponseData.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseData.weights/max(thePhotocurrentImpulseResponseData.weights), 'ko-');
    hold on;
    p2 = plot(theHighPassImpulseResponseData.temporalSupportSeconds*1e3, theHighPassImpulseResponseData.weights/max(theHighPassImpulseResponseData.weights), 'b--');
    p3 = plot(theMacaqueCenterImpulseResponseData.temporalSupportSeconds*1e3, highPassPrediction/max(highPassPrediction), 'ms-');
    legend([p1, p2, p3], {'pCurrent', 'high-pass', 'high-pass prediction'})
    xlabel('time (msec)')
    set(gca, 'XLim', [0 500], 'XTick', 0:10:1000);
    grid on


    pause
    % Fit
    gain = struct(...
        'initial', centerParams(1), ...
        'low', centerParams(1)/100, ...
        'high', centerParams(1)*100 ...
        );

    highPassGain = struct(...
        'initial',centerParams(2), ...
        'low', 0, ...
        'high', 10 ...
        );

    highPassTimeConstant = struct(...
        'initial',centerParams(3), ...
        'low', 0, ...
        'high', 100 ...
        );

    lowPassStagesNum = struct(...
        'initial',centerParams(4), ...
        'low', 1, ...
        'high', 300 ...
        );

    lowPassTimeConstant = struct(...
        'initial',centerParams(5), ...
        'low', 0.1, ...
        'high', 200 ...
        );

    delay = struct(...
        'initial', centerParams(6), ...
        'low', 0, ...
        'high', 400 ...
        );

    TTFparams.initialValues = [gain.initial   highPassGain.initial    highPassTimeConstant.initial    lowPassStagesNum.initial  lowPassTimeConstant.initial      delay.initial];
    TTFparams.lowerBounds   = [gain.low       highPassGain.low        highPassTimeConstant.low        lowPassStagesNum.low      lowPassTimeConstant.low          delay.low];
    TTFparams.upperBounds   = [gain.high      highPassGain.high       highPassTimeConstant.high       lowPassStagesNum.high     lowPassTimeConstant.high         delay.high];
    TTFparams.names         = {'gain',        'highPassGain',         'high pass tau (msec)',  'low pass stages',       'low pass tau (msec)',    'delay (msec)'};
    TTFparams.scaling       = {'log',         'linear',               'linear',                       'linear',                 'linear',                        'linear'};

    hFig = figure(100);
    axModelParams = subplot(1,1,1);

    % Fit theProductTTF
    [TTFparams, theFitttedTTF] = highPassNstageLowPassTTFtoArbitraryTTF(...
        theComplexTTF, temporalFrequenciesExamined, TTFparams, axModelParams);

    % Check the fit
    figure(33); clf;
    axTTFamplitude = subplot(1,1,1);
    p1 = plot(axTTFamplitude, temporalFrequenciesExamined, theTTFamplitude, 'ko-', 'MarkerSize', 12);
    hold on;
    p2 = plot(axTTFamplitude, temporalFrequenciesExamined, abs(theFitttedTTF), 'r.-', 'LineWidth', 1.5);
    legend(axTTFamplitude, [p1 p2], {'data', 'fit'});
    set(axTTFamplitude, 'XLim', [0.3 300], 'XTick', [0.3 1 3 10 30 100 300], 'XScale', 'log');
    grid(axTTFamplitude, 'on');
    xlabel(axTTFamplitude,'frequency (Hz)');
    ylabel(axTTFamplitude,'amplitude');


    % Convert to temporal impulse response functions

    theFittedImpulseResponseData = ...
        temporalTransferFunctionToImpulseResponseFunction(theFitttedTTF, temporalFrequenciesExamined);

    % Plot impulse response functions
    figure(60); clf;
    plot(theFittedImpulseResponseData.temporalSupportSeconds*1e3, theFittedImpulseResponseData.weights, 'k.-');
    xlabel('time (msec)')
    set(gca, 'XTick', 0:10:1000);
    grid on

end





function theImpulseResponseFunctionData = ...
        temporalTransferFunctionToImpulseResponseFunction(theTTF, temporalFrequencySupportHz, ...
        zeroPaddingLength, delayMilliseconds, performFFTshift)


    % Add delay in frequency domain
    theTTF  = exp(-1i*2*pi*temporalFrequencySupportHz*delayMilliseconds/1000) .* theTTF;

    % Zero padding
    if (~isempty(zeroPaddingLength))
        extraSamplesNum = 512-numel(temporalFrequencySupportHz);
        nSamples = numel(theTTF) + extraSamplesNum;
        theTTF = cat(2, theTTF, zeros(1,nSamples-numel(theTTF)));
        temporalFrequencySupportHz = temporalFrequencySupportHz(1) + (0:(nSamples-1))*(temporalFrequencySupportHz(2)-temporalFrequencySupportHz(1));
    end

    % Convert single-sided spectrum to double sided
    theDoubleSidedTTF = [theTTF(1) theTTF(2:end)/2 fliplr(conj(theTTF(2:end)))/2];

    if (performFFTshift)
        theImpulseResponse = fftshift(ifft(theDoubleSidedTTF, 'symmetric'));
    else
        theImpulseResponse = ifft(theDoubleSidedTTF, 'symmetric');
    end

    fMax = max(temporalFrequencySupportHz);
    dtSeconds = 1/(2*fMax)
    theTemporalSupportSeconds = (1:numel(theImpulseResponse)) * dtSeconds;

    theImpulseResponseFunctionData.weights = theImpulseResponse;
    theImpulseResponseFunctionData.temporalSupportSeconds = theTemporalSupportSeconds;
end


% Fit the highPassNstageLowPassTTF to an arbitrary TTF
function  [TTFparams, theFittedTTF] = highPassNstageLowPassTTFtoArbitraryTTF(theArbitraryComplexTTF, temporalFrequencySupportHz, TTFparams, ax)

    objective = @(x)highPassNstageLowPassTTFresidual(x, theArbitraryComplexTTF, temporalFrequencySupportHz, ax, TTFparams);

    % Multi-start
    problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', TTFparams.initialValues, ...
          'lb', TTFparams.lowerBounds, ...
          'ub', TTFparams.upperBounds, ...
          'options', optimoptions(...
            'fmincon',...
            'Display', 'none', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^4) ...
          );

    ms = MultiStart(...
          'Display', 'iter', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', ~true);

    multiStartsNum = 128;

    % Run the multi-start
    TTFparams.finalValues = run(ms, problem, multiStartsNum);


    theFittedTTF = highPassNstageLowPassTTF(TTFparams.finalValues, temporalFrequencySupportHz);

end


function theTTF = highPassNstageTTF(params, temporalFrequencySupportHz)

    % Get params
    gain = params(1);
    highPassStagesNum = round(params(2));
    highPassTimeConstantSeconds = params(3)*1e-3;
    delaySeconds = params(4)*1e-3;
    omega = temporalFrequencySupportHz * (2 * pi);

    % The TwoStageTTF model in the frequency domain
    j_omega_tau = 1i * omega * highPassTimeConstantSeconds;
    theTTF = ...
        gain * exp(-1i * omega * delaySeconds) .* ((j_omega_tau .^ highPassStagesNum) ./ ((1 + j_omega_tau) .^ highPassStagesNum));


end

function theTTF = highPassNstageLowPassTTF(params, temporalFrequencySupportHz)

    % Get params
    gain = params(1);
    highPassGain = params(2);
    highPassTimeConstantSeconds = params(3)*1e-3;
    lowPassStagesNum = round(params(4));
    lowPassTimeConstantSeconds = params(5)*1e-3;
    delaySeconds = params(6)*1e-3;
    omega = temporalFrequencySupportHz * (2 * pi);

    % The TwoStageTTF model in the frequency domain
    theTTF = ...
        gain * exp(-1i * omega * delaySeconds) .* ...
        (1 - highPassGain * 1./(1 + 1i * omega * highPassTimeConstantSeconds)) .* ...
        (1./(1 + 1i * omega * lowPassTimeConstantSeconds)) .^ (lowPassStagesNum);

end


function theResidual = highPassNstageLowPassTTFresidual(theCurrentParams, theTTFtoFit, temporalFrequencySupportHz, ax, modelVariables)

    theResidual = norm(highPassNstageLowPassTTF(theCurrentParams, temporalFrequencySupportHz) - theTTFtoFit);

    modelVariables.finalValues = theCurrentParams;
    RGCMosaicConstructor.visualize.fittedModelParams(ax, modelVariables, 'TTF fit');

end



