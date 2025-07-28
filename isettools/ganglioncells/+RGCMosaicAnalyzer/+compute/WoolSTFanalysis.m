function WoolSTFanalysis(...
    theMRGCMosaicSTFResponsesFullFileName, ...
    theWoolAnalysisFileName, ...
    theLowSF, ...
    targetedSurroundPurityRange, ...
    targetedRadialEccentricityRange, ...
    targetedCenterConeNumerosityRange, ...
    targetedCenterPurityRange, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('fixedOptimalOrientation', [], @(x)(ischar(x)||isnumeric(x)||(isnan(x))));
    p.addParameter('visualizeSingleRGCConeIsolatingSTFs', false, @islogical);
    p.addParameter('singleRGCindexToAnalyze', [], @(x)(isscalar(x)||isempty(x)));
    p.addParameter('classifyCenterDominanceBasedOnLowFrequencyResponse', false, @islogical);
    p.addParameter('classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights', true, @islogical);
    p.addParameter('incrementExcitatoryLowSFresponsePhaseDegs', 270, @isscalar);
    p.addParameter('computeSurroundConePurityExcludingCenterConnectedCones', false, @islogical);
    p.addParameter('unwrapSTFphase', true, @islogical);
    p.addParameter('unwrapLMopponentSTFphase', true, @islogical);
    p.addParameter('forceBandPassSTFphaseUnwrap', ~true, @islogical);
    p.addParameter('employNativeSTFphaseUnwrapMethod', true, @islogical);
    p.addParameter('maxResponseStrengthRatio', 1.3, @isscalar);
    p.addParameter('reAnalyzeData', true, @islogical);

    p.parse(varargin{:});
    classifyCenterDominanceBasedOnLowFrequencyResponse = p.Results.classifyCenterDominanceBasedOnLowFrequencyResponse;
    classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights = p.Results.classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights;
    fixedOptimalOrientation = p.Results.fixedOptimalOrientation;
    visualizeSingleRGCConeIsolatingSTFs = p.Results.visualizeSingleRGCConeIsolatingSTFs;
    singleRGCindexToAnalyze = p.Results.singleRGCindexToAnalyze;
    incrementExcitatoryLowSFresponsePhaseDegs = p.Results.incrementExcitatoryLowSFresponsePhaseDegs;
    computeSurroundConePurityExcludingCenterConnectedCones = p.Results.computeSurroundConePurityExcludingCenterConnectedCones;
    employNativeSTFphaseUnwrapMethod = p.Results.employNativeSTFphaseUnwrapMethod;
    forceBandPassSTFphaseUnwrap = p.Results.forceBandPassSTFphaseUnwrap;
    unwrapSTFphase = p.Results.unwrapSTFphase;
    unwrapLMopponentSTFphase = p.Results.unwrapLMopponentSTFphase;
    maxResponseStrengthRatio = p.Results.maxResponseStrengthRatio;
    reAnalyzeData = p.Results.reAnalyzeData;

    % Where PDF files are exported
    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'WoolResponseStrengthLMphase');
    thePDFfileName2 = fullfile(theRawFiguresDir, pdfExportSubDir, 'WoolConePurity');

    % Bin cell eccentricities in bins of 0.2 mm
    binnedEccMMs = 0:0.2:8;

    % Allow up to 32 cells in each of these bins
    maxCellsPerBin = 32;

    if (~reAnalyzeData)
        generateAllRunPlots(binnedEccMMs, maxCellsPerBin, ...
            maxResponseStrengthRatio, thePDFfileName, thePDFfileName2);
        return;
    end

    if (classifyCenterDominanceBasedOnLowFrequencyResponse && classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights)
        error('You cannot set both ''classifyCenterDominanceBasedOnLowFrequencyResponse'' and ''classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights'' to true.\nOne of them must be set to true, and the other one to false');
    end
    if (~classifyCenterDominanceBasedOnLowFrequencyResponse && ~classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights)
        error('You cannot set both ''classifyCenterDominanceBasedOnLowFrequencyResponse'' and ''classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights'' to false.\nOne of them must be set to true, and the other one to false');
    end


    if (iscell(theWoolAnalysisFileName))
        runsNum = numel(theWoolAnalysisFileName);
    else
        runsNum = 1;
    end

    responseStrengthRatiosLcenterDominatedAllRuns = [];
    responseStrengthRatiosMcenterDominatedAllRuns = [];
    phaseDifferenceDegsLcenterDominatedAllRuns = [];
    phaseDifferenceDegsMcenterDominatedAllRuns = [];
    

    centerConePurityOpponentLcentersAllRuns = [];
    centerConePurityOpponentMcentersAllRuns =  [];
    surroundConePurityOpponentLcentersAllRuns =  [];
    surroundConePurityOpponentMcentersAllRuns =  [];
        

    centerConePurityNonOpponentLcentersAllRuns =  [];
    centerConePurityNonOpponentMcentersAllRuns =  [];
    surroundConePurityNonOpponentLcentersAllRuns =  [];
    surroundConePurityNonOpponentMcentersAllRuns =  [];
        

    opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = [];
    opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = [];
    nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = [];
    nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = [];

    opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = [];
    opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = [];
    nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = [];
    nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = [];

    opponentLcenterDominatedEccMMsAllRuns = [];
    opponentMcenterDominatedEccMMsAllRuns = [];
    nonopponentLcenterDominatedEccMMsAllRuns = [];
    nonopponentMcenterDominatedEccMMsAllRuns = [];



    runHorizontalEccMMs = zeros(1, runsNum);

    for iRun = 1:runsNum

        if (iscell(theWoolAnalysisFileName))
            theCurrentWoolAnalysisFileName = theWoolAnalysisFileName{iRun};
            theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName{iRun};
        else
            theCurrentWoolAnalysisFileName = theWoolAnalysisFileName;
            theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName;
        end

        if (isempty(strfind(theCurrentMRGCMosaicSTFResponsesFullFileName, 'LconeIsolating')))
            theMRGCMosaicMconeIsolatingSTFResponsesFullFileName = theCurrentMRGCMosaicSTFResponsesFullFileName;
            theMRGCMosaicLconeIsolatingSTFResponsesFullFileName = strrep(theMRGCMosaicMconeIsolatingSTFResponsesFullFileName, ...
                'MconeIsolating', 'LconeIsolating');
        else
            theMRGCMosaicLconeIsolatingSTFResponsesFullFileName = theCurrentMRGCMosaicSTFResponsesFullFileName;
            theMRGCMosaicMconeIsolatingSTFResponsesFullFileName = strrep(theMRGCMosaicLconeIsolatingSTFResponsesFullFileName, ...
                'LconeIsolating', 'MconeIsolating');
        end

        % Load the L-cone isolating responses and the mosaic
		load(theMRGCMosaicLconeIsolatingSTFResponsesFullFileName, ...
            'theMRGCMosaic', ...
            'stimParams', ...
			'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF', 'theMRGCMosaicResponseTemporalSupportSeconds');
        LconeIsolatingSpatioTemporalMRGCMosaicResponses2DSTF = theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF;

        % Load the M-cone isolating responses
        load(theMRGCMosaicMconeIsolatingSTFResponsesFullFileName, ...
			'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF');
        MconeIsolatingSpatioTemporalMRGCMosaicResponses2DSTF = theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF;

        LcenterDominatedLowSFConeIsolatingSTFdata = [];
        McenterDominatedLowSFConeIsolatingSTFdata = [];

        LcenterDominatedPeakSFConeIsolatingSTFdata = [];
        McenterDominatedPeakSFConeIsolatingSTFdata = [];

        responseStrengthRatiosLcenterDominated = [];
        responseStrengthRatiosMcenterDominated = [];
        phaseDifferenceDegsLcenterDominated = [];
        phaseDifferenceDegsMcenterDominated = [];

        LcenterDominatedIndicesBasedOnLowSFresponse = [];
        McenterDominatedIndicesBasedOnLowSFresponse = [];

        LcenterDominatedIndicesBasedOnRFcenterPoolingWeights = [];
        McenterDominatedIndicesBasedOnRFcenterPoolingWeights = [];
        
        opponentLcenterIdxBasedOnRFpoolingWeights = [];
        opponentMcenterIdxBasedOnRFpoolingWeights = [];
        nonopponentLcenterIdxBasedOnRFpoolingWeights = [];
        nonopponentMcenterIdxBasedOnRFpoolingWeights = [];

        centerConePurity = zeros(1, theMRGCMosaic.rgcsNum);
        surroundConePurity = zeros(1, theMRGCMosaic.rgcsNum);

        for iRGC = 1:theMRGCMosaic.rgcsNum

            if (~isempty(singleRGCindexToAnalyze))
                if (iRGC ~= singleRGCindexToAnalyze)
                    continue;
                end
            end

            fprintf('Analyzing RGC %d of %d\n', iRGC, theMRGCMosaic.rgcsNum);

            % Get the analyzed L-- and M--cone isolating STF slices for this cell
            [theLconeIsolatingSTF, theMconeIsolatingSTF, sfSupport, centerConeDominance] = extractTheAnalyzedSTFslices(...
                theMRGCMosaic, iRGC, stimParams, ...
                theMRGCMosaicResponseTemporalSupportSeconds, ...
                squeeze(LconeIsolatingSpatioTemporalMRGCMosaicResponses2DSTF(:,:,:,iRGC)), ...
				squeeze(MconeIsolatingSpatioTemporalMRGCMosaicResponses2DSTF(:,:,:,iRGC)), ...
                fixedOptimalOrientation, ...
                incrementExcitatoryLowSFresponsePhaseDegs);

            % What frequencies to use as 'low SF'
            if (numel(theLowSF) == 1)
                [~,lowFrequencyIndices] = min(abs(sfSupport - theLowSF));
            elseif (numel(theLowSF) == 2) && (theLowSF(2)>=theLowSF(1))
                lowFrequencyIndices = find((sfSupport >= theLowSF(1)) & (sfSupport <= theLowSF(2)));
            else
                error('theLowSF must have either 1 (single SF) or 2 elements (range of SFs)')
            end

            % Compute the non-dominant cone to dominant cone response
            % strength ratio vs the L-M response phase difference using the
            % low spatial frequencies
            lowFrequencyLconeResponseAmplitude = mean(theLconeIsolatingSTF.amplitudeSpectrum(lowFrequencyIndices));
            lowFrequencyMconeResponseAmplitude = mean(theMconeIsolatingSTF.amplitudeSpectrum(lowFrequencyIndices));
            lowFrequencyLconeResponsePhase     = mean(theLconeIsolatingSTF.phaseSpectrum(lowFrequencyIndices));
            lowFrequencyMconeResponsePhase     = mean(theMconeIsolatingSTF.phaseSpectrum(lowFrequencyIndices));


            % Assign signs to the amplitudes. 
            % %Phase near 0/360 degs (+-45 degs)
            % are assigned a positive sign (response excitatory to increment).
            % Phases near +- 180 (+/- 45 degs) 
            % are assigned a negative sign (response excitatory to decrement).
            if (abs(lowFrequencyLconeResponsePhase-0) < 45) || (abs(lowFrequencyLconeResponsePhase-360) < 45)
                lowFrequencyLconeResponseAmplitudeSign = 1;
            elseif (abs(abs(lowFrequencyLconeResponsePhase)-180) < 45)
                lowFrequencyLconeResponseAmplitudeSign = -1;
            else
                if (abs(lowFrequencyLconeResponsePhase-0) < 90) || (abs(lowFrequencyLconeResponsePhase-360) < 90)
                    lowFrequencyLconeResponseAmplitudeSign = 1;
                else
                    lowFrequencyLconeResponseAmplitudeSign = -1;
                end
                fprintf(2,'L-cone response phase is not within +/- 45 from 0/360 or 180/-180 (%f)\n.', lowFrequencyLconeResponsePhase);
            end

            % Check the response at the phase that corresponds to the increment part of the M-cone stimulus
            if (abs(lowFrequencyMconeResponsePhase-0) < 45) || (abs(lowFrequencyMconeResponsePhase-360) < 45)
                lowFrequencyMconeResponseAmplitudeSign = 1;
            elseif (abs(abs(lowFrequencyMconeResponsePhase)-180) < 45)
                lowFrequencyMconeResponseAmplitudeSign = -1;
            else
                if (abs(lowFrequencyMconeResponsePhase-0) <= 90) || (abs(lowFrequencyMconeResponsePhase-360) <= 90)
                    lowFrequencyMconeResponseAmplitudeSign = 1;
                else
                    lowFrequencyMconeResponseAmplitudeSign = -1;
                end
                fprintf(2,'M-cone response phase is not within +/- 45 0/360 or 180/-180(%f)\n.', lowFrequencyMconeResponsePhase);
            end

            theLastLcenterDominatedResponseStrengthRatio = [];
            theLastLcenterDominatedPhaseDifference = [];
            theLastMcenterDominatedResponseStrengthRatio = [];
            theLastMcenterDominatedPhaseDifference = [];

            if (lowFrequencyLconeResponseAmplitude*lowFrequencyLconeResponseAmplitudeSign > lowFrequencyMconeResponseAmplitude*lowFrequencyMconeResponseAmplitudeSign) 
                % Response strength ratio: non-dominant cone type /
                % dominant cone type at the lowest SF
                theLastLcenterDominatedResponseStrengthRatio = lowFrequencyMconeResponseAmplitude/lowFrequencyLconeResponseAmplitude;
                theLastLcenterDominatedPhaseDifference = lowFrequencyLconeResponsePhase - lowFrequencyMconeResponsePhase;

                LcenterDominatedIndicesBasedOnLowSFresponse(numel(LcenterDominatedIndicesBasedOnLowSFresponse)+1) = iRGC;
                responseStrengthRatiosLcenterDominated(numel(responseStrengthRatiosLcenterDominated)+1) = theLastLcenterDominatedResponseStrengthRatio;
                phaseDifferenceDegsLcenterDominated(numel(phaseDifferenceDegsLcenterDominated)+1) = theLastLcenterDominatedPhaseDifference;

                LcenterDominatedLowSFConeIsolatingSTFdata(size(LcenterDominatedLowSFConeIsolatingSTFdata,1)+1,:) = [ ...
                    lowFrequencyLconeResponseAmplitude ...
                    lowFrequencyMconeResponseAmplitude ...
                    lowFrequencyLconeResponsePhase ...
                    lowFrequencyMconeResponsePhase ];

                % Determine peak spatial from the L-cone STF
                [~,peakSFindex] = max(theLconeIsolatingSTF.amplitudeSpectrum(:));

            else
                theLastMcenterDominatedResponseStrengthRatio = lowFrequencyLconeResponseAmplitude/lowFrequencyMconeResponseAmplitude;
                theLastMcenterDominatedPhaseDifference = lowFrequencyLconeResponsePhase - lowFrequencyMconeResponsePhase;

                McenterDominatedIndicesBasedOnLowSFresponse(numel(McenterDominatedIndicesBasedOnLowSFresponse)+1) = iRGC;
                responseStrengthRatiosMcenterDominated(numel(responseStrengthRatiosMcenterDominated)+1)= theLastMcenterDominatedResponseStrengthRatio;
                phaseDifferenceDegsMcenterDominated(numel(phaseDifferenceDegsMcenterDominated)+1) = theLastMcenterDominatedPhaseDifference;
                    
                McenterDominatedLowSFConeIsolatingSTFdata(size(McenterDominatedLowSFConeIsolatingSTFdata,1)+1,:) = [ ...
                    lowFrequencyLconeResponseAmplitude ...
                    lowFrequencyMconeResponseAmplitude ...
                    lowFrequencyLconeResponsePhase ...
                    lowFrequencyMconeResponsePhase ];

                % Determine peak spatial frequency from the M-cone STF
                [~,peakSFindex] = max(theMconeIsolatingSTF.amplitudeSpectrum(:));
            end

            if (lowFrequencyLconeResponseAmplitude*lowFrequencyLconeResponseAmplitudeSign > lowFrequencyMconeResponseAmplitude*lowFrequencyMconeResponseAmplitudeSign)
                LcenterDominatedPeakSFConeIsolatingSTFdata(size(LcenterDominatedPeakSFConeIsolatingSTFdata,1)+1,:) = [ ...
                    theLconeIsolatingSTF.amplitudeSpectrum(peakSFindex) ...
                    theMconeIsolatingSTF.amplitudeSpectrum(peakSFindex) ...
                    theLconeIsolatingSTF.phaseSpectrum(peakSFindex) ...
                    theMconeIsolatingSTF.phaseSpectrum(peakSFindex)];
            else
                McenterDominatedPeakSFConeIsolatingSTFdata(size(McenterDominatedPeakSFConeIsolatingSTFdata,1)+1,:) = [ ...
                    theLconeIsolatingSTF.amplitudeSpectrum(peakSFindex) ...
                    theMconeIsolatingSTF.amplitudeSpectrum(peakSFindex) ...
                    theLconeIsolatingSTF.phaseSpectrum(peakSFindex) ...
                    theMconeIsolatingSTF.phaseSpectrum(peakSFindex)];
            end

            % Wool et al compute cL, cM, sL, sM (and therefore cone purity) from fitting the DoG model to the
            % L-- and M-- cone isolating STFs at 2 Hz. This probably produces the variation they see.
            % We dont have a temporal component in
            % our mRGC model, so we compute cL, cM, sL, sM as the spatially
            % integrated cone weights in the center and surround subregions
            % We could also fit the measured STF but for very low spatial
            % frequencies and no optics, this would be equivalent to the
            % spatiailly integrated L and M cone weights, which is what we
            % do here.

            % Compute cone and surround L/(L+M) purity indices
            theCenterConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
                iRGC, 'center', ...
                'minConeWeightIncluded', 0.0001, ...
                'inputConeIndicesOnly', true, ...
                'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
            theCenterConeWeights = full(squeeze(theMRGCMosaic.rgcRFcenterConeConnectivityMatrix(theCenterConeIndices, iRGC)));

            if (computeSurroundConePurityExcludingCenterConnectedCones)
                theSurroundConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
                    iRGC, 'surround-center', ...
                    'minConeWeightIncluded', 0.0001, ...
                    'inputConeIndicesOnly', true, ...
                    'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
            else
                theSurroundConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
                    iRGC, 'surround', ...
                    'minConeWeightIncluded', 0.0001, ...
                    'inputConeIndicesOnly', true, ...
                    'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);
            end
            theSurroundConeWeights = full(squeeze(theMRGCMosaic.rgcRFsurroundConeConnectivityMatrix(theSurroundConeIndices, iRGC)));

            % Spatially integrated L cone weight in the RF center
            idxL = find(theMRGCMosaic.inputConeMosaic.coneTypes(theCenterConeIndices) == cMosaic.LCONE_ID);
            cL = sum(theCenterConeWeights(idxL));

            % Spatially integrated M cone weight in the RF center
            idxM = find(theMRGCMosaic.inputConeMosaic.coneTypes(theCenterConeIndices) == cMosaic.MCONE_ID);
            cM = sum(theCenterConeWeights(idxM));

            % Spatially integrated L cone weight in the RF surround
            idxL = find(theMRGCMosaic.inputConeMosaic.coneTypes(theSurroundConeIndices) == cMosaic.LCONE_ID);
            sL = sum(theSurroundConeWeights(idxL));

            % Spatially integrated M cone weight in the RF surround
            idxM = find(theMRGCMosaic.inputConeMosaic.coneTypes(theSurroundConeIndices) == cMosaic.MCONE_ID);
            sM = sum(theSurroundConeWeights(idxM));

            % Compute cone purities from spatially integrated cone weights
            % to center and surrounds
            centerConePurity(iRGC) = cL/(cL+cM);
            surroundConePurity(iRGC) = sL/(sL+sM);


            if (classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights)
                % Spatially integrated L-cone weight over the entire RF
                totalLweight = cL-sL;
    
                % Spatially integrated M-cone weight over the entire RF
                totalMweight = cM-sM;

                if (cL > cM)
                    theIndex = numel(LcenterDominatedIndicesBasedOnRFcenterPoolingWeights)+1;
                    LcenterDominatedIndicesBasedOnRFcenterPoolingWeights(theIndex) = iRGC;
                   
                    if (totalLweight * totalMweight < 0)
                        % Opponent L-center dominated cell
                        opponentLcenterIdxBasedOnRFpoolingWeights(numel(opponentLcenterIdxBasedOnRFpoolingWeights)+1) = theIndex;
                    else
                        % Non-opponent L-center dominated cell
                        nonopponentLcenterIdxBasedOnRFpoolingWeights(numel(nonopponentLcenterIdxBasedOnRFpoolingWeights)+1) = theIndex;
                    end
    
                else
                    theIndex = numel(McenterDominatedIndicesBasedOnRFcenterPoolingWeights)+1;
                    McenterDominatedIndicesBasedOnRFcenterPoolingWeights(theIndex) = iRGC;
                    
                    if (totalLweight * totalMweight < 0)
                        % Opponent M-center dominated cell
                        opponentMcenterIdxBasedOnRFpoolingWeights(numel(opponentMcenterIdxBasedOnRFpoolingWeights)+1) = theIndex;
                    else
                        % Non-opponent L-center dominated cell
                        nonopponentMcenterIdxBasedOnRFpoolingWeights(numel(nonopponentMcenterIdxBasedOnRFpoolingWeights)+1) = theIndex;
                    end
                end
            end %  if (classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights)

            

            if (visualizeSingleRGCConeIsolatingSTFs)

                % Visualize STF and RF
                hFig = visualizeSingleCellConeWeightsMapsAndConeIsolatingSTFs(...
                    theMRGCMosaic, iRGC, centerConeDominance, sfSupport, ...
                    theLconeIsolatingSTF, theMconeIsolatingSTF, ...
                    unwrapSTFphase, employNativeSTFphaseUnwrapMethod, forceBandPassSTFphaseUnwrap, ...
                    unwrapLMopponentSTFphase); 

                % Export plot
                theSingleCellConeIsolatingSTFsPDFfileName = ...
                    fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('singleCellConeIsolatingSTF_%2.2f_%2.2f_RGC%d.pdf', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, iRGC));
                NicePlot.exportFigToPDF(theSingleCellConeIsolatingSTFsPDFfileName, hFig,  300);


                % Generate the responseRatioPhaseDifferent plot for this cell
                [hFig, opponentLcenterIndex, nonOpponentLcenterIndex, ...
                       opponentMcenterIndex, nonOpponentMcenterIndex] = generateResponseRatioPhasePlot(...
                    50, ...
                    theLastLcenterDominatedResponseStrengthRatio, ...
                    theLastMcenterDominatedResponseStrengthRatio, ...
                    theLastLcenterDominatedPhaseDifference, ...
                    theLastMcenterDominatedPhaseDifference, ...
                    maxResponseStrengthRatio, ...
                    sprintf('ecc, xy (mm): (%2.2f, %2.2f), RGC: %d/%d', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, iRGC, theMRGCMosaic.rgcsNum));

                % Export it
                theSingleCellResponseRatioPhasePlotPDFfileName = ...
                    fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('singleCellResponseRatioPhasePlot_%2.2f_%2.2f_RGC%d.pdf', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, iRGC));
                NicePlot.exportFigToPDF(theSingleCellResponseRatioPhasePlotPDFfileName, hFig,  300);

                % Generate the cone purity plot for this cell
                if (~isempty(opponentLcenterIndex))
                    hFig = generateConePurityPlot(...
                        100, ...
                        centerConePurity(iRGC), surroundConePurity(iRGC), ...
                        [], [], ...
                        [], [], ...
                        [], [], ...
                        'all', ...
                        sprintf('ecc, xy (mm): (%2.2f, %2.2f), RGC: %d/%d', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0,iRGC, theMRGCMosaic.rgcsNum));
                        
                elseif (~isempty(opponentMcenterIndex))
                    hFig = generateConePurityPlot(...
                        101, ...
                        [], [], ...
                        centerConePurity(iRGC), surroundConePurity(iRGC), ...
                        [], [], ...
                        [], [], ...
                        'all', ...
                        sprintf('ecc, xy (mm): (%2.2f, %2.2f), RGC: %d/%d', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0,iRGC, theMRGCMosaic.rgcsNum));
                        

                elseif (~isempty(nonOpponentLcenterIndex))
                    hFig = generateConePurityPlot(...
                        102, ...
                        [], [], ...
                        [], [], ...
                        centerConePurity(iRGC), surroundConePurity(iRGC), ...
                        [], [], ...
                        'all', ...
                        sprintf('ecc, xy (mm): (%2.2f, %2.2f), RGC: %d/%d', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0,iRGC, theMRGCMosaic.rgcsNum));
                
                elseif (~isempty(nonOpponentMcenterIndex))
                    hFig = generateConePurityPlot(...
                        103, ...
                        [], [], ...
                        [], [], ...
                        [], [], ...
                        centerConePurity(iRGC), surroundConePurity(iRGC), ...
                        'all', ...
                        sprintf('ecc, xy (mm): (%2.2f, %2.2f), RGC: %d/%d', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0,iRGC, theMRGCMosaic.rgcsNum));
                end
 
                % Export it
                theSingleCellResponseRatioPhasePlotPDFfileName = ...
                    fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('singleCellConePurityPlot_%2.2f_%2.2f_RGC%d.pdf', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, iRGC));
                NicePlot.exportFigToPDF(theSingleCellResponseRatioPhasePlotPDFfileName, hFig,  300);


            end

        end % iRGC



        if (1==2)
            % Generate the LowSF response/phase plot
            hFig = generateConeIsolatingLowSFresponsedataPlot(999, ...
                LcenterDominatedLowSFConeIsolatingSTFdata, ...
                McenterDominatedLowSFConeIsolatingSTFdata);
    
            pause
        end

        % Generate the responseRatioPhaseDifferent plot for this run
        [hFig,opponentLcenterIdxBasedOnLowSFresponse, nonopponentLcenterIdxBasedOnLowSFresponse, ...
              opponentMcenterIdxBasedOnLowSFresponse, nonopponentMcenterIdxBasedOnLowSFresponse] = generateResponseRatioPhasePlot(...
                500, ...
                responseStrengthRatiosLcenterDominated, ...
                responseStrengthRatiosMcenterDominated, ...
                phaseDifferenceDegsLcenterDominated, ...
                phaseDifferenceDegsMcenterDominated, ...
                maxResponseStrengthRatio, ...
                sprintf('ecc, xy (mm): (%2.2f, %2.2f)', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0));

        % Export it
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end

        NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


        LcenterDominatedEccMMs = 1e-3 * sqrt(sum(theMRGCMosaic.rgcRFpositionsMicrons(LcenterDominatedIndicesBasedOnLowSFresponse,:).^2,2));
        McenterDominatedEccMMs = 1e-3 * sqrt(sum(theMRGCMosaic.rgcRFpositionsMicrons(McenterDominatedIndicesBasedOnLowSFresponse,:).^2,2));

        % Split the data into L- and M-center dominated opponent cells and L- and M-center dominated non-opponent cells
        if (classifyCenterDominanceBasedOnLowFrequencyResponse)
            opponentLcenterIndices = LcenterDominatedIndicesBasedOnLowSFresponse(opponentLcenterIdxBasedOnLowSFresponse);
            opponentMcenterIndices = McenterDominatedIndicesBasedOnLowSFresponse(opponentMcenterIdxBasedOnLowSFresponse);
            nonOpponentLcenterIndices = LcenterDominatedIndicesBasedOnLowSFresponse(nonopponentLcenterIdxBasedOnLowSFresponse);
            nonOpponentMcenterIndices = McenterDominatedIndicesBasedOnLowSFresponse(nonopponentMcenterIdxBasedOnLowSFresponse);

            opponentLcenterDominatedLowSFConeIsolatingSTFdata = LcenterDominatedLowSFConeIsolatingSTFdata(opponentLcenterIdxBasedOnLowSFresponse,:);
            opponentMcenterDominatedLowSFConeIsolatingSTFdata = McenterDominatedLowSFConeIsolatingSTFdata(opponentMcenterIdxBasedOnLowSFresponse,:);
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata = LcenterDominatedLowSFConeIsolatingSTFdata(nonopponentLcenterIdxBasedOnLowSFresponse,:);
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata = McenterDominatedLowSFConeIsolatingSTFdata(nonopponentMcenterIdxBasedOnLowSFresponse,:);

            opponentLcenterDominatedEccMMs = LcenterDominatedEccMMs(opponentLcenterIdxBasedOnLowSFresponse);
            opponentMcenterDominatedEccMMs = McenterDominatedEccMMs(opponentMcenterIdxBasedOnLowSFresponse);
            nonOpponentLcenterDominatedEccMMs = LcenterDominatedEccMMs(nonopponentLcenterIdxBasedOnLowSFresponse);
            nonOpponentMcenterDominatedEccMMs = McenterDominatedEccMMs(nonopponentMcenterIdxBasedOnLowSFresponse);

            % Do the same for LcenterDominatedPeakSFConeIsolatingSTFdata
            opponentLcenterDominatedPeakSFConeIsolatingSTFdata = LcenterDominatedPeakSFConeIsolatingSTFdata(opponentLcenterIdxBasedOnLowSFresponse,:);
            opponentMcenterDominatedPeakSFConeIsolatingSTFdata = McenterDominatedPeakSFConeIsolatingSTFdata(opponentMcenterIdxBasedOnLowSFresponse,:);
            nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata = LcenterDominatedPeakSFConeIsolatingSTFdata(nonopponentLcenterIdxBasedOnLowSFresponse,:);
            nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata = McenterDominatedPeakSFConeIsolatingSTFdata(nonopponentMcenterIdxBasedOnLowSFresponse,:);


        elseif (classifyCenterDominanceBasedOnRFcenterIntegratedPoolingWeights)
            opponentLcenterIndices = LcenterDominatedIndicesBasedOnRFcenterPoolingWeights(opponentLcenterIdxBasedOnRFpoolingWeights);
            opponentMcenterIndices = McenterDominatedIndicesBasedOnRFcenterPoolingWeights(opponentMcenterIdxBasedOnRFpoolingWeights);
            nonOpponentLcenterIndices = LcenterDominatedIndicesBasedOnRFcenterPoolingWeights(nonopponentLcenterIdxBasedOnRFpoolingWeights);
            nonOpponentMcenterIndices = McenterDominatedIndicesBasedOnRFcenterPoolingWeights(nonopponentMcenterIdxBasedOnRFpoolingWeights);

            opponentLcenterDominatedLowSFConeIsolatingSTFdata = LcenterDominatedLowSFConeIsolatingSTFdata(opponentLcenterIdxBasedOnRFpoolingWeights,:);
            opponentMcenterDominatedLowSFConeIsolatingSTFdata = McenterDominatedLowSFConeIsolatingSTFdata(opponentMcenterIdxBasedOnRFpoolingWeights,:);
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata = LcenterDominatedLowSFConeIsolatingSTFdata(nonopponentLcenterIdxBasedOnRFpoolingWeights,:);
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata = McenterDominatedLowSFConeIsolatingSTFdata(nonopponentMcenterIdxBasedOnRFpoolingWeights,:);

            opponentLcenterDominatedEccMMs = LcenterDominatedEccMMs(opponentLcenterIdxBasedOnRFpoolingWeights);
            opponentMcenterDominatedEccMMs = McenterDominatedEccMMs(opponentMcenterIdxBasedOnRFpoolingWeights);
            nonOpponentLcenterDominatedEccMMs = LcenterDominatedEccMMs(nonopponentLcenterIdxBasedOnRFpoolingWeights);
            nonOpponentMcenterDominatedEccMMs = McenterDominatedEccMMs(nonopponentMcenterIdxBasedOnRFpoolingWeights);

            % Do the same for LcenterDominatedPeakSFConeIsolatingSTFdata
            opponentLcenterDominatedPeakSFConeIsolatingSTFdata = LcenterDominatedPeakSFConeIsolatingSTFdata(opponentLcenterIdxBasedOnRFpoolingWeights,:);
            opponentMcenterDominatedPeakSFConeIsolatingSTFdata = McenterDominatedPeakSFConeIsolatingSTFdata(opponentMcenterIdxBasedOnRFpoolingWeights,:);
            nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata = LcenterDominatedPeakSFConeIsolatingSTFdata(nonopponentLcenterIdxBasedOnRFpoolingWeights,:);
            nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata = McenterDominatedPeakSFConeIsolatingSTFdata(nonopponentMcenterIdxBasedOnRFpoolingWeights,:);

        end


        % Generate plot of lowSF L-M opponent response as a function of cell's eccentricity in MMs
        hFig = generateLowSFLMopponentResponseVariationWithEccentricity(888, ...
            opponentLcenterDominatedLowSFConeIsolatingSTFdata, ...
            opponentMcenterDominatedLowSFConeIsolatingSTFdata, ...
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata, ...
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata, ...
            opponentLcenterDominatedEccMMs, ...
            opponentMcenterDominatedEccMMs, ...
            nonOpponentLcenterDominatedEccMMs, ...
            nonOpponentMcenterDominatedEccMMs);

        % Export it
        runHorizontalEccMMs(iRun) = theMRGCMosaic.eccentricityMicrons(1)/1000.0;
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_LMopponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_LMopponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


        % Generate plot of peak SF L+M non-opponent response as a function of cell's eccentricity in MMs
        hFig = generatePeakSFLMnonOpponentResponseVariationWithEccentricity(889, ...
             opponentLcenterDominatedPeakSFConeIsolatingSTFdata, ...
             opponentMcenterDominatedPeakSFConeIsolatingSTFdata, ...
             nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata, ...
             nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata, ...
             opponentLcenterDominatedEccMMs, ...
             opponentMcenterDominatedEccMMs, ...
             nonOpponentLcenterDominatedEccMMs, ...
             nonOpponentMcenterDominatedEccMMs);

        % Export it
        runHorizontalEccMMs(iRun) = theMRGCMosaic.eccentricityMicrons(1)/1000.0;
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_LMnonOpponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_LMnonOpponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


        % Generate the cone purity plot for all cells this run
        [hFig,hFigCenterPurityHistogram, hFigSurroundPurityHistogram] = generateConePurityPlot(...
            200, ...
            centerConePurity(opponentLcenterIndices), surroundConePurity(opponentLcenterIndices), ...
            centerConePurity(opponentMcenterIndices), surroundConePurity(opponentMcenterIndices), ...
            centerConePurity(nonOpponentLcenterIndices), surroundConePurity(nonOpponentLcenterIndices), ...
            centerConePurity(nonOpponentMcenterIndices), surroundConePurity(nonOpponentMcenterIndices), ...
            'opponentOnly', ...
            sprintf('ecc, xy (mm): (%2.2f, %2.2f)', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0) ...
            );

        % Export the scatter
        runHorizontalEccMMs(iRun) = theMRGCMosaic.eccentricityMicrons(1)/1000.0;
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_Opponent_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_Opponent_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


        % Export the center histogram
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_OpponentCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_OpponentCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFigCenterPurityHistogram,  300);

        % Export the surround histogram
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_OpponentSurroundHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_OpponentSurroundCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFigSurroundPurityHistogram,  300);


        [hFig,hFigCenterPurityHistogram, hFigSurroundPurityHistogram] = generateConePurityPlot(...
            300, ...
            centerConePurity(opponentLcenterIndices), surroundConePurity(opponentLcenterIndices), ...
            centerConePurity(opponentMcenterIndices), surroundConePurity(opponentMcenterIndices), ...
            centerConePurity(nonOpponentLcenterIndices), surroundConePurity(nonOpponentLcenterIndices), ...
            centerConePurity(nonOpponentMcenterIndices), surroundConePurity(nonOpponentMcenterIndices), ...
            'nonOpponentOnly', ...
            sprintf('ecc, xy (mm): (%2.2f, %2.2f)', theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0) ...
            );

        % Export the scatter
        runHorizontalEccMMs(iRun) = theMRGCMosaic.eccentricityMicrons(1)/1000.0;
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_NonOpponent_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_NonOpponent_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);

        % Export the center histogram
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_NonOpponentCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_NonOpponentCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFigCenterPurityHistogram,  300);

        % Export the surround histogram
        if (numel(theLowSF) == 1)
            theRunPDFfileName = sprintf('%s_NonOpponentSurroundHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF);
        else
            theRunPDFfileName = sprintf('%s_NonOpponentSurroundCenterHistogram_EccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, theMRGCMosaic.eccentricityMicrons(1)/1000.0, theMRGCMosaic.eccentricityMicrons(2)/1000.0, theLowSF(1), theLowSF(2));
        end
        NicePlot.exportFigToPDF(theRunPDFfileName,hFigSurroundPurityHistogram,  300);


        % Concatenate data from this run to the all runs
        centerConePurityOpponentLcentersAllRuns = cat(1, centerConePurityOpponentLcentersAllRuns(:), centerConePurity(opponentLcenterIndices)');
        centerConePurityOpponentMcentersAllRuns = cat(1, centerConePurityOpponentMcentersAllRuns(:), centerConePurity(opponentMcenterIndices)');
        surroundConePurityOpponentLcentersAllRuns = cat(1, surroundConePurityOpponentLcentersAllRuns(:), surroundConePurity(opponentLcenterIndices)');
        surroundConePurityOpponentMcentersAllRuns = cat(1, surroundConePurityOpponentMcentersAllRuns(:), surroundConePurity(opponentMcenterIndices)');
        
        centerConePurityNonOpponentLcentersAllRuns = cat(1, centerConePurityNonOpponentLcentersAllRuns(:), centerConePurity(nonOpponentLcenterIndices)');
        centerConePurityNonOpponentMcentersAllRuns = cat(1, centerConePurityNonOpponentMcentersAllRuns(:), centerConePurity(nonOpponentMcenterIndices)');
        surroundConePurityNonOpponentLcentersAllRuns = cat(1, surroundConePurityNonOpponentLcentersAllRuns(:), surroundConePurity(nonOpponentLcenterIndices)');
        surroundConePurityNonOpponentMcentersAllRuns = cat(1, surroundConePurityNonOpponentMcentersAllRuns(:), surroundConePurity(nonOpponentMcenterIndices)');
        
        responseStrengthRatiosLcenterDominatedAllRuns = cat(1, responseStrengthRatiosLcenterDominatedAllRuns(:), responseStrengthRatiosLcenterDominated(:));
        responseStrengthRatiosMcenterDominatedAllRuns = cat(1, responseStrengthRatiosMcenterDominatedAllRuns(:), responseStrengthRatiosMcenterDominated(:));
        phaseDifferenceDegsLcenterDominatedAllRuns = cat(1, phaseDifferenceDegsLcenterDominatedAllRuns(:), phaseDifferenceDegsLcenterDominated(:));
        phaseDifferenceDegsMcenterDominatedAllRuns = cat(1, phaseDifferenceDegsMcenterDominatedAllRuns(:), phaseDifferenceDegsMcenterDominated(:));



        opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = cat(1, ...
            opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns, opponentLcenterDominatedLowSFConeIsolatingSTFdata);
        opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = cat(1, ...
            opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns, opponentMcenterDominatedLowSFConeIsolatingSTFdata);
        nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = cat(1, ...
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns, nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata);
        nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = cat(1, ...
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns, nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata);


        opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = cat(1, ...
            opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, opponentLcenterDominatedPeakSFConeIsolatingSTFdata);
        opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = cat(1, ...
            opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, opponentMcenterDominatedPeakSFConeIsolatingSTFdata);
        nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = cat(1, ...
            nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata);
        nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = cat(1, ...
            nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata);

        opponentLcenterDominatedEccMMsAllRuns = cat(1, opponentLcenterDominatedEccMMsAllRuns(:), opponentLcenterDominatedEccMMs);
        opponentMcenterDominatedEccMMsAllRuns = cat(1, opponentMcenterDominatedEccMMsAllRuns(:), opponentMcenterDominatedEccMMs);
        nonopponentLcenterDominatedEccMMsAllRuns = cat(1, nonopponentLcenterDominatedEccMMsAllRuns(:), nonOpponentLcenterDominatedEccMMs);
        nonopponentMcenterDominatedEccMMsAllRuns = cat(1, nonopponentMcenterDominatedEccMMsAllRuns(:), nonOpponentMcenterDominatedEccMMs);
    end % iRun

    % Save analyzed data from all runs
    theAllRunsAnalysisFileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'WoolAnalysisAllRuns.mat');

    save(theAllRunsAnalysisFileName, ...
        'responseStrengthRatiosLcenterDominatedAllRuns', ...
        'responseStrengthRatiosMcenterDominatedAllRuns', ...
        'phaseDifferenceDegsLcenterDominatedAllRuns', ...
        'phaseDifferenceDegsMcenterDominatedAllRuns', ...
        'centerConePurityOpponentLcentersAllRuns', 'surroundConePurityOpponentLcentersAllRuns', ...
        'centerConePurityOpponentMcentersAllRuns', 'surroundConePurityOpponentMcentersAllRuns', ...
        'centerConePurityNonOpponentLcentersAllRuns', 'surroundConePurityNonOpponentLcentersAllRuns', ...
        'centerConePurityNonOpponentMcentersAllRuns', 'surroundConePurityNonOpponentMcentersAllRuns', ...
        'opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'opponentLcenterDominatedEccMMsAllRuns', ...
        'opponentMcenterDominatedEccMMsAllRuns', ...
        'nonopponentLcenterDominatedEccMMsAllRuns', ...
        'nonopponentMcenterDominatedEccMMsAllRuns', ...
        'runHorizontalEccMMs', 'theLowSF', ...
        '-v7.3');

    generateAllRunPlots(binnedEccMMs, maxCellsPerBin, ...
            maxResponseStrengthRatio, thePDFfileName, thePDFfileName2);
end


function generateAllRunPlots(binnedEccMMs, maxCellsPerBin, maxResponseStrengthRatio, thePDFfileName, thePDFfileName2)

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    theAllRunsAnalysisFileName = fullfile(theRawFiguresDir, pdfExportSubDir, 'WoolAnalysisAllRuns.mat');

    load(theAllRunsAnalysisFileName, ...
        'responseStrengthRatiosLcenterDominatedAllRuns', ...
        'responseStrengthRatiosMcenterDominatedAllRuns', ...
        'phaseDifferenceDegsLcenterDominatedAllRuns', ...
        'phaseDifferenceDegsMcenterDominatedAllRuns', ...
        'centerConePurityOpponentLcentersAllRuns', 'surroundConePurityOpponentLcentersAllRuns', ...
        'centerConePurityOpponentMcentersAllRuns', 'surroundConePurityOpponentMcentersAllRuns', ...
        'centerConePurityNonOpponentLcentersAllRuns', 'surroundConePurityNonOpponentLcentersAllRuns', ...
        'centerConePurityNonOpponentMcentersAllRuns', 'surroundConePurityNonOpponentMcentersAllRuns', ...
        'opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns', ...
        'opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns', ...
        'opponentLcenterDominatedEccMMsAllRuns', ...
        'opponentMcenterDominatedEccMMsAllRuns', ...
        'nonopponentLcenterDominatedEccMMsAllRuns', ...
        'nonopponentMcenterDominatedEccMMsAllRuns', ...
        'runHorizontalEccMMs', 'theLowSF');

    % Subsample the opponent L-center dominated cells
    idx = subsampleCells(100, 'opponent, L-center', opponentLcenterDominatedEccMMsAllRuns, binnedEccMMs, maxCellsPerBin);
    opponentLcenterDominatedEccMMsAllRuns = opponentLcenterDominatedEccMMsAllRuns(idx);
    centerConePurityOpponentLcentersAllRuns = centerConePurityOpponentLcentersAllRuns(idx);
    surroundConePurityOpponentLcentersAllRuns = surroundConePurityOpponentLcentersAllRuns(idx);
    opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns(idx,:);
    opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns(idx,:)

    % Subsample the opponent M-center dominated cells
    idx = subsampleCells(101, 'opponent, M-center',opponentMcenterDominatedEccMMsAllRuns, binnedEccMMs, maxCellsPerBin);
    opponentMcenterDominatedEccMMsAllRuns = opponentMcenterDominatedEccMMsAllRuns(idx);
    centerConePurityOpponentMcentersAllRuns = centerConePurityOpponentMcentersAllRuns(idx);
    surroundConePurityOpponentMcentersAllRuns = surroundConePurityOpponentMcentersAllRuns(idx);
    opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns(idx,:);
    opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns(idx,:);


    % Subsample the non-opponent L-center dominated cells
    idx = subsampleCells(102, 'non-opponent, L-center', nonopponentLcenterDominatedEccMMsAllRuns, binnedEccMMs, maxCellsPerBin);
    nonopponentLcenterDominatedEccMMsAllRuns = nonopponentLcenterDominatedEccMMsAllRuns(idx);
    centerConePurityNonOpponentLcentersAllRuns = centerConePurityNonOpponentLcentersAllRuns(idx);
    surroundConePurityNonOpponentLcentersAllRuns = surroundConePurityNonOpponentLcentersAllRuns(idx);
    nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns = nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns(idx,:);
    nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns(idx,:);

    % Subsample the non-opponent M-center dominated cells
    idx = subsampleCells(103, 'non-opponent, M-center', nonopponentMcenterDominatedEccMMsAllRuns, binnedEccMMs, maxCellsPerBin);
    nonopponentMcenterDominatedEccMMsAllRuns = nonopponentMcenterDominatedEccMMsAllRuns(idx);
    centerConePurityNonOpponentMcentersAllRuns = centerConePurityNonOpponentMcentersAllRuns(idx);
    surroundConePurityNonOpponentMcentersAllRuns = surroundConePurityNonOpponentMcentersAllRuns(idx);
    nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns = nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns(idx,:);
    nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns = nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns(idx,:);


    % Generate the responseRatioPhaseDifferent plot across all analyzed runs
    hFig = generateResponseRatioPhasePlot(...
                5000,...
                responseStrengthRatiosLcenterDominatedAllRuns, ...
                responseStrengthRatiosMcenterDominatedAllRuns, ...
                phaseDifferenceDegsLcenterDominatedAllRuns, ...
                phaseDifferenceDegsMcenterDominatedAllRuns, ...
                maxResponseStrengthRatio, ...
                sprintf('ecc, x (mm): %2.2f - %2.2f', min(runHorizontalEccMMs), max(runHorizontalEccMMs)));

    % Export it
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


    % Generate the cone purity plot for opponent cells across  all analyzed runs
    [hFig,hFigCenterPurityHistogram, hFigSurroundPurityHistogram] = generateConePurityPlot(...
            1000, ...
            centerConePurityOpponentLcentersAllRuns, surroundConePurityOpponentLcentersAllRuns, ...
            centerConePurityOpponentMcentersAllRuns, surroundConePurityOpponentMcentersAllRuns, ...
            centerConePurityNonOpponentLcentersAllRuns, surroundConePurityNonOpponentLcentersAllRuns, ...
            centerConePurityNonOpponentMcentersAllRuns, surroundConePurityNonOpponentMcentersAllRuns, ...
            'opponentOnly', ...
            sprintf('ecc, x (mm): %2.2f - %2.2f', min(runHorizontalEccMMs), max(runHorizontalEccMMs)));

    % Export the scatter
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_Opponent_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_Opponent_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


    % Export the center histogram
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_OpponentCenterHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_OpponentCenterHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFigCenterPurityHistogram,  300);


    % Export the surround histogram
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_OpponentSurroundHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_OpponentSurroundHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFigSurroundPurityHistogram,  300);



    % Generate the cone purity plot for non-opponent cells across  all analyzed runs
    [hFig,hFigCenterPurityHistogram, hFigSurroundPurityHistogram] = generateConePurityPlot(...
            2000, ...
            centerConePurityOpponentLcentersAllRuns, surroundConePurityOpponentLcentersAllRuns, ...
            centerConePurityOpponentMcentersAllRuns, surroundConePurityOpponentMcentersAllRuns, ...
            centerConePurityNonOpponentLcentersAllRuns, surroundConePurityNonOpponentLcentersAllRuns, ...
            centerConePurityNonOpponentMcentersAllRuns, surroundConePurityNonOpponentMcentersAllRuns, ...
            'nonOpponentOnly', ...
            sprintf('ecc, x (mm): %2.2f - %2.2f', min(runHorizontalEccMMs), max(runHorizontalEccMMs)));

    % Export it
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_NonOpponent_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_NonOpponent_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);


    % Export the center histogram
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_NonOpponentCenterHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_NonOpponentCenterHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFigCenterPurityHistogram,  300);


    % Export the surround histogram
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_NonOpponentSurroundHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_NonOpponentSurroundHistogram_EccMMs_%2.2f_to_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
                thePDFfileName2, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFigSurroundPurityHistogram,  300);


    % Generate plot of lowSF L-M opponent response as a function of cell's eccentricity in MMs
    hFig = generateLowSFLMopponentResponseVariationWithEccentricity(888, ...
            opponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns, ...
            opponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns, ...
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdataAllRuns, ...
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdataAllRuns, ...
            opponentLcenterDominatedEccMMsAllRuns, ...
            opponentMcenterDominatedEccMMsAllRuns, ...
            nonopponentLcenterDominatedEccMMsAllRuns, ...
            nonopponentMcenterDominatedEccMMsAllRuns);


    % Export it
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_LMopponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
            thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_LMopponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
            thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end
    NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);



    % Generate plot of peak SF L+M non-opponent response as a function of cell's eccentricity in MMs
    hFig = generatePeakSFLMnonOpponentResponseVariationWithEccentricity(889, ...
            opponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, ...
            opponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, ...
            nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, ...
            nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdataAllRuns, ...
            opponentLcenterDominatedEccMMsAllRuns, ...
            opponentMcenterDominatedEccMMsAllRuns, ...
            nonopponentLcenterDominatedEccMMsAllRuns, ...
            nonopponentMcenterDominatedEccMMsAllRuns);


    % Export it
    if (numel(theLowSF) == 1)
        theRunPDFfileName = sprintf('%s_LMnonOpponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2fCPD.pdf', ...
            thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF);
    else
        theRunPDFfileName = sprintf('%s_LMnonPpponentResponseVsEccMMs_%2.2f_%2.2f_LowSF_%2.2f-%2.2fCPD.pdf', ...
            thePDFfileName, min(runHorizontalEccMMs), max(runHorizontalEccMMs), theLowSF(1), theLowSF(2));
    end

    NicePlot.exportFigToPDF(theRunPDFfileName,hFig,  300);

end


function cellIndicesToKeep = subsampleCells(figNo, cellType, eccMMsAllRuns, binnedEccMMs, maxCellsPerBin)

    [N,edges,bins] = histcounts(eccMMsAllRuns, binnedEccMMs);

    cellsAtBin = zeros(1, max(bins));
    cellIndicesAtBin = cell(1, max(bins));
    for iBin = 1:max(bins)
        % How many cells in this ecc bin
        theCellIndicesAtThisEccBin = find(bins == iBin);
        cellsAtBin(iBin) = numel(theCellIndicesAtThisEccBin);
        cellIndicesAtBin{iBin} = theCellIndicesAtThisEccBin;
    end

    cellIndicesToKeep = [];
    for iBin = 1:max(bins)
        newCellIndices = cellIndicesAtBin{iBin};
        if (cellsAtBin(iBin)>maxCellsPerBin)
            % Only keep maxCellsPerBin
            idx = randperm(numel(newCellIndices));
            newCellIndices = newCellIndices(idx(1:maxCellsPerBin));
        end
        fprintf('Cells at bin %d: %d\n', iBin, numel(newCellIndices));
        cellIndicesToKeep = cat(1, cellIndicesToKeep(:), newCellIndices(:));
    end

    figure(figNo); clf;
    ax = subplot(1,2,1)
    histogram(ax, eccMMsAllRuns, binnedEccMMs);
    xlabel(ax, 'ecc (mm)');
    ylabel(ax, 'cell count');
    title(ax, sprintf('ecc distribution all %s cells', cellType));

    ax = subplot(1,2,2)
    histogram(ax, eccMMsAllRuns(cellIndicesToKeep), binnedEccMMs);
    xlabel(ax, 'ecc (mm)');
    ylabel(ax, 'cell count');
    title(ax, sprintf('ecc distribution for subsampled %s cells', cellType));
    drawnow;
end


function hFig = generatePeakSFLMnonOpponentResponseVariationWithEccentricity(figNo, ...
            opponentLcenterDominatedPeakSFConeIsolatingSTFdata, ...
            opponentMcenterDominatedPeakSFConeIsolatingSTFdata, ...
            nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata, ...
            nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata, ...
            opponentLcenterDominatedEccMMs, ...
            opponentMcenterDominatedEccMMs, ...
            nonOpponentLcenterDominatedEccMMs, ...
            nonOpponentMcenterDominatedEccMMs)

    % The L+M non opponent responses of opponent L-center dominated RGCs
    opponentLcenterDominatedLconeSTFamplitude = opponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,1);
    opponentLcenterDominatedMconeSTFamplitude = opponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,2);
    opponentLcenterDominatedLconeSTFphase     = opponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,3);
    opponentLcenterDominatedMconeSTFphase     = opponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,4);

    opponentLcenterDominatedLMnonOpponentResponses = abs(...
        opponentLcenterDominatedLconeSTFamplitude .* exp(1i*opponentLcenterDominatedLconeSTFphase/ 180*pi) + ...
        opponentLcenterDominatedMconeSTFamplitude .* exp(1i*opponentLcenterDominatedMconeSTFphase / 180*pi));


     % The L+M non opponent responses of opponent M-center dominated RGCs
    opponentMcenterDominatedLconeSTFamplitude = opponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,1);
    opponentMcenterDominatedMconeSTFamplitude = opponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,2);
    opponentMcenterDominatedLconeSTFphase     = opponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,3);
    opponentMcenterDominatedMconeSTFphase     = opponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,4);

    opponentMcenterDominatedLMnonOpponentResponses = abs(...
        opponentMcenterDominatedLconeSTFamplitude .* exp(1i*opponentMcenterDominatedLconeSTFphase/ 180*pi) + ...
        opponentMcenterDominatedMconeSTFamplitude .* exp(1i*opponentMcenterDominatedMconeSTFphase / 180*pi));

    % The L+M non-opponent responses of non-opponent L-center dominated RGCs
    nonOpponentLcenterDominatedLconeSTFamplitude = nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,1);
    nonOpponentLcenterDominatedMconeSTFamplitude = nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,2);
    nonOpponentLcenterDominatedLconeSTFphase     = nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,3);
    nonOpponentLcenterDominatedMconeSTFphase     = nonOpponentLcenterDominatedPeakSFConeIsolatingSTFdata(:,4);

    nonOpponentLcenterDominatedLMnonOpponentResponses = abs(...
        nonOpponentLcenterDominatedLconeSTFamplitude .* exp(1i*nonOpponentLcenterDominatedLconeSTFphase/ 180*pi) + ...
        nonOpponentLcenterDominatedMconeSTFamplitude .* exp(1i*nonOpponentLcenterDominatedMconeSTFphase / 180*pi));

    % The L+M non-opponent responses of non-opponent M-center dominated RGCs
    nonOpponentMcenterDominatedLconeSTFamplitude = nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,1);
    nonOpponentMcenterDominatedMconeSTFamplitude = nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,2);
    nonOpponentMcenterDominatedLconeSTFphase     = nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,3);
    nonOpponentMcenterDominatedMconeSTFphase     = nonOpponentMcenterDominatedPeakSFConeIsolatingSTFdata(:,4);

    nonOpponentMcenterDominatedLMnonOpponentResponses = abs(...
        nonOpponentMcenterDominatedLconeSTFamplitude .* exp(1i*nonOpponentMcenterDominatedLconeSTFphase/ 180*pi) + ...
        nonOpponentMcenterDominatedMconeSTFamplitude .* exp(1i*nonOpponentMcenterDominatedMconeSTFphase / 180*pi));

    maxResponse = max([...
        max(opponentLcenterDominatedLMnonOpponentResponses) ...
        max(opponentMcenterDominatedLMnonOpponentResponses) ...
        max(nonOpponentLcenterDominatedLMnonOpponentResponses) ...
        max(nonOpponentMcenterDominatedLMnonOpponentResponses)]);
    maxResponse = 1.2;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    markerSize = ff.markerSize-8;
    lineWidth = ff.lineWidth/2;

    scatter(ax, opponentLcenterDominatedEccMMs, opponentLcenterDominatedLMnonOpponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor',  LconeSTFcolor*0.5, 'MarkerFaceColor',  LconeSTFcolor);

    hold(ax, 'on')
    scatter(ax, opponentMcenterDominatedEccMMs, opponentMcenterDominatedLMnonOpponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor',  MconeSTFcolor*0.5, 'MarkerFaceColor',  MconeSTFcolor);

    scatter(ax, nonOpponentLcenterDominatedEccMMs, nonOpponentLcenterDominatedLMnonOpponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
         'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);

    scatter(ax, nonOpponentMcenterDominatedEccMMs, nonOpponentMcenterDominatedLMnonOpponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
         'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);

    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [0 7], 'YLim', [0 maxResponse], 'XTick', 0:1:10, 'YTick', 0:0.2:2);

    xlabel(ax, 'eccentricity (mm)');
    ylabel(ax, 'L+M non-opponent response');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    drawnow;

end



function hFig = generateLowSFLMopponentResponseVariationWithEccentricity(figNo, ...
            opponentLcenterDominatedLowSFConeIsolatingSTFdata, ...
            opponentMcenterDominatedLowSFConeIsolatingSTFdata, ...
            nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata, ...
            nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata, ...
            opponentLcenterDominatedEccMMs, ...
            opponentMcenterDominatedEccMMs, ...
            nonOpponentLcenterDominatedEccMMs, ...
            nonOpponentMcenterDominatedEccMMs)

    % The L-M opponent responses of opponent L-center dominated RGCs
    opponentLcenterDominatedLconeSTFamplitude = opponentLcenterDominatedLowSFConeIsolatingSTFdata(:,1);
    opponentLcenterDominatedMconeSTFamplitude = opponentLcenterDominatedLowSFConeIsolatingSTFdata(:,2);
    opponentLcenterDominatedLconeSTFphase     = opponentLcenterDominatedLowSFConeIsolatingSTFdata(:,3);
    opponentLcenterDominatedMconeSTFphase     = opponentLcenterDominatedLowSFConeIsolatingSTFdata(:,4);

    opponentLcenterDominatedLMopponentResponses = abs(...
        opponentLcenterDominatedLconeSTFamplitude .* exp(1i*opponentLcenterDominatedLconeSTFphase/ 180*pi) - ...
        opponentLcenterDominatedMconeSTFamplitude .* exp(1i*opponentLcenterDominatedMconeSTFphase / 180*pi));


    % The L-M opponent responses of opponent M-center dominated RGCs
    opponentMcenterDominatedLconeSTFamplitude = opponentMcenterDominatedLowSFConeIsolatingSTFdata(:,1);
    opponentMcenterDominatedMconeSTFamplitude = opponentMcenterDominatedLowSFConeIsolatingSTFdata(:,2);
    opponentMcenterDominatedLconeSTFphase     = opponentMcenterDominatedLowSFConeIsolatingSTFdata(:,3);
    opponentMcenterDominatedMconeSTFphase     = opponentMcenterDominatedLowSFConeIsolatingSTFdata(:,4);

    opponentMcenterDominatedLMopponentResponses = abs(...
        opponentMcenterDominatedLconeSTFamplitude .* exp(1i*opponentMcenterDominatedLconeSTFphase/ 180*pi) - ...
        opponentMcenterDominatedMconeSTFamplitude .* exp(1i*opponentMcenterDominatedMconeSTFphase / 180*pi));


    % The L-M opponent responses of non-opponent L-center dominated RGCs
    nonOpponentLcenterDominatedLconeSTFamplitude = nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata(:,1);
    nonOpponentLcenterDominatedMconeSTFamplitude = nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata(:,2);
    nonOpponentLcenterDominatedLconeSTFphase     = nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata(:,3);
    nonOpponentLcenterDominatedMconeSTFphase     = nonOpponentLcenterDominatedLowSFConeIsolatingSTFdata(:,4);

    nonOpponentLcenterDominatedLMopponentResponses = abs(...
        nonOpponentLcenterDominatedLconeSTFamplitude .* exp(1i*nonOpponentLcenterDominatedLconeSTFphase/ 180*pi) - ...
        nonOpponentLcenterDominatedMconeSTFamplitude .* exp(1i*nonOpponentLcenterDominatedMconeSTFphase / 180*pi));

    % The L-M opponent responses of non-opponent M-center dominated RGCs
    nonOpponentMcenterDominatedLconeSTFamplitude = nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata(:,1);
    nonOpponentMcenterDominatedMconeSTFamplitude = nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata(:,2);
    nonOpponentMcenterDominatedLconeSTFphase     = nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata(:,3);
    nonOpponentMcenterDominatedMconeSTFphase     = nonOpponentMcenterDominatedLowSFConeIsolatingSTFdata(:,4);

    nonOpponentMcenterDominatedLMopponentResponses = abs(...
        nonOpponentMcenterDominatedLconeSTFamplitude .* exp(1i*nonOpponentMcenterDominatedLconeSTFphase/ 180*pi) - ...
        nonOpponentMcenterDominatedMconeSTFamplitude .* exp(1i*nonOpponentMcenterDominatedMconeSTFphase / 180*pi));



    maxResponse = max([...
        max(opponentLcenterDominatedLMopponentResponses) ...
        max(opponentMcenterDominatedLMopponentResponses) ...
        max(nonOpponentLcenterDominatedLMopponentResponses) ...
        max(nonOpponentMcenterDominatedLMopponentResponses)]);
    maxResponse = 1.2;

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    markerSize = ff.markerSize-8;
    lineWidth = ff.lineWidth/2;

    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    scatter(ax, opponentLcenterDominatedEccMMs, opponentLcenterDominatedLMopponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor', LconeSTFcolor*0.5, 'MarkerFaceColor', LconeSTFcolor);

    hold(ax, 'on')
    scatter(ax, opponentMcenterDominatedEccMMs, opponentMcenterDominatedLMopponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor', MconeSTFcolor*0.5, 'MarkerFaceColor', MconeSTFcolor);

    scatter(ax, nonOpponentLcenterDominatedEccMMs, nonOpponentLcenterDominatedLMopponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
         'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);

    scatter(ax, nonOpponentMcenterDominatedEccMMs, nonOpponentMcenterDominatedLMopponentResponses, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
         'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);

    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [0 7], 'YLim', [0 maxResponse], 'XTick', 0:1:10, 'YTick', 0:0.2:2);

    xlabel(ax, 'eccentricity (mm)');
    ylabel(ax, 'L-M opponent response');

    PublicationReadyPlotLib.applyFormat(ax,ff);

    drawnow;

end


function hFig = generateConeIsolatingLowSFresponsedataPlot(figNo,...
    LcenterDominatedLowSFConeIsolatingSTFdata, McenterDominatedLowSFConeIsolatingSTFdata)

    lowFrequencyLconeResponseAmplitudesLdom = LcenterDominatedLowSFConeIsolatingSTFdata(:,1);
    lowFrequencyMconeResponseAmplitudesLdom = LcenterDominatedLowSFConeIsolatingSTFdata(:,2);
    lowFrequencyLconeResponsePhasesLdom     = LcenterDominatedLowSFConeIsolatingSTFdata(:,3);
    lowFrequencyMconeResponsePhasesLdom     = LcenterDominatedLowSFConeIsolatingSTFdata(:,4);

    lowFrequencyLconeResponseAmplitudesMdom = McenterDominatedLowSFConeIsolatingSTFdata(:,1);
    lowFrequencyMconeResponseAmplitudesMdom = McenterDominatedLowSFConeIsolatingSTFdata(:,2);
    lowFrequencyLconeResponsePhasesMdom     = McenterDominatedLowSFConeIsolatingSTFdata(:,3);
    lowFrequencyMconeResponsePhasesMdom     = McenterDominatedLowSFConeIsolatingSTFdata(:,4);


    figure(444);
    ax = subplot(2,2,1);
    phases = [lowFrequencyLconeResponsePhasesLdom; lowFrequencyMconeResponsePhasesLdom];
    amplitudes = [lowFrequencyLconeResponseAmplitudesLdom; lowFrequencyMconeResponseAmplitudesLdom];
    plot(phases, amplitudes, 'k-');
    hold (ax, 'on')
    p1 = plot(ax,lowFrequencyLconeResponsePhasesLdom, lowFrequencyLconeResponseAmplitudesLdom, 'ro');
    p2 = plot(ax,lowFrequencyMconeResponsePhasesLdom, lowFrequencyMconeResponseAmplitudesLdom, 'go');
    set(ax, 'XLim', [0 360], 'YLim', [0 1]);
    legend(ax, [p1 p2], {'L-cone isolating responses', 'M-cone isolating responses'})
    xlabel(ax, 'response phase (degs)');
    ylabel(ax, 'response amplitude');
    title(ax, 'L-center dominated');

    ax = subplot(2,2,2);
    phases = [lowFrequencyLconeResponsePhasesMdom; lowFrequencyMconeResponsePhasesMdom];
    amplitudes = [lowFrequencyLconeResponseAmplitudesMdom; lowFrequencyMconeResponseAmplitudesMdom];
    plot(phases, amplitudes, 'k-');
    hold (ax, 'on')
    p1 = plot(ax,lowFrequencyLconeResponsePhasesMdom, lowFrequencyLconeResponseAmplitudesMdom, 'ro');
    p2 = plot(ax,lowFrequencyMconeResponsePhasesMdom, lowFrequencyMconeResponseAmplitudesMdom, 'go');
    set(ax, 'XLim', [0 360], 'YLim', [0 1]);
    legend(ax, [p1 p2], {'L-cone isolating responses', 'M-cone isolating responses'})
    xlabel(ax, 'response phase (degs)');
    ylabel(ax, 'response amplitude');
    title(ax, 'M-center dominated');


    pause


    maxResponse = max([...
        max(lowFrequencyLconeResponseAmplitudesLdom(:)) ...
        max(lowFrequencyMconeResponseAmplitudesLdom(:)) ...
        max(lowFrequencyLconeResponseAmplitudesMdom(:)) ...
        max(lowFrequencyMconeResponseAmplitudesMdom(:)) ...
        ]);


    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    % The L-center dominated cells

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    markerSize = ff.markerSize-8;
    lineWidth = ff.lineWidth/2;
    
    % Render the polar grid
    theRadii = (0.2:0.2:1.0)*maxResponse;
    theAngles = 0:30:360;
    renderPolarGrid(ax, theRadii, theAngles, ff.lineWidth/2, [0.4 0.4 0.4], 'L-center dominated cells');

    xL = lowFrequencyLconeResponseAmplitudesLdom .* cosd(lowFrequencyLconeResponsePhasesLdom);
    yL = lowFrequencyLconeResponseAmplitudesLdom .* sind(lowFrequencyLconeResponsePhasesLdom);
    xM = lowFrequencyMconeResponseAmplitudesLdom .* cosd(lowFrequencyMconeResponsePhasesLdom);
    yM = lowFrequencyMconeResponseAmplitudesLdom .* sind(lowFrequencyMconeResponsePhasesLdom);

    
    scatter(ax, xL,yL, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
        'MarkerEdgeColor', LconeSTFcolor*0.5, 'MarkerFaceColor', LconeSTFcolor);

    hold(ax, 'on')

    scatter(ax, xM,yM, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
        'MarkerEdgeColor', LconeSTFcolor*0.5, 'MarkerFaceColor', LconeSTFcolor);

    axis(ax,'equal')
    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [-1 1]*maxResponse, 'YLim', [-1 1]*maxResponse, 'XTick', [], 'YTick', []);
    
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'XColor', 'none', 'YColor', 'none');

    drawnow;


    % The M-center dominated cells
    hFig = figure(figNo+1); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    % Render the polar grid
    theRadii = (0.2:0.2:1.0)*maxResponse;
    theAngles = 0:30:360;
    renderPolarGrid(ax, theRadii, theAngles, maxResponse, ff.lineWidth/2, [0.4 0.4 0.4], 'M-center dominated cells');

    xL = lowFrequencyLconeResponseAmplitudesMdom .* cosd(lowFrequencyLconeResponsePhasesMdom);
    yL = lowFrequencyLconeResponseAmplitudesMdom .* sind(lowFrequencyLconeResponsePhasesMdom);
    xM = lowFrequencyMconeResponseAmplitudesMdom .* cosd(lowFrequencyMconeResponsePhasesMdom);
    yM = lowFrequencyMconeResponseAmplitudesMdom .* sind(lowFrequencyMconeResponsePhasesMdom);

    
    scatter(ax, xL,yL, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
        'MarkerEdgeColor', MconeSTFcolor*0.5, 'MarkerFaceColor', MconeSTFcolor);

    hold(ax, 'on')

    scatter(ax, xM,yM, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
        'MarkerEdgeColor', MconeSTFcolor*0.5, 'MarkerFaceColor', MconeSTFcolor);


    axis(ax,'equal')
    set(ax, 'FontSize', ff.axisFontSize);
    
    set(ax, 'XLim', [-1 1]*maxResponse, 'YLim', [-1 1]*maxResponse, 'XTick', [], 'YTick', []);
    
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'XColor', 'none', 'YColor', 'none');

    drawnow;
end

function [hFigPurityScatter, hFigCenterPurityHistogram, hFigSurroundPurityHistogram] = generateConePurityPlot( ...
             figNo, ...
             centerConePurityOpponentLcenters, surroundConePurityOpponentLcenters, ...
             centerConePurityOpponentMcenters, surroundConePurityOpponentMcenters, ...
             centerConePurityNonOpponentLcenters, surroundConePurityNonOpponentLcenters, ...
             centerConePurityNonOpponentMcenters, surroundConePurityNonOpponentMcenters, ...
             whichCells, ...
             plotTitle)


    
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFigPurityScatter = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFigPurityScatter,ff);
    ax = theAxes{1,1};

    if (max([...
        numel(centerConePurityOpponentLcenters)
        numel(centerConePurityOpponentMcenters)
        numel(centerConePurityNonOpponentLcenters)
        numel(centerConePurityNonOpponentMcenters)]) == 1)
        markerSize = ff.markerSize;
        lineWidth = ff.lineWidth;
    else
        markerSize = ff.markerSize-8;
        lineWidth = ff.lineWidth/2;
    end

    hold(ax, 'on')
    
    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    if (~isempty(centerConePurityOpponentLcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'opponentOnly')))
        scatter(ax, centerConePurityOpponentLcenters, surroundConePurityOpponentLcenters, markerSize^2, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
            'MarkerEdgeColor', LconeSTFcolor*0.5, 'MarkerFaceColor', LconeSTFcolor);
    end
    
    if (~isempty(centerConePurityOpponentMcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'opponentOnly')))
        scatter(ax, centerConePurityOpponentMcenters, surroundConePurityOpponentMcenters, markerSize^2, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
            'MarkerEdgeColor',  MconeSTFcolor*0.5, 'MarkerFaceColor',MconeSTFcolor);
    end

    if (~isempty(centerConePurityNonOpponentLcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'nonOpponentOnly')))
        scatter(ax, centerConePurityNonOpponentLcenters, surroundConePurityNonOpponentLcenters, markerSize^2, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.8 0.8 0.8]);
    end

    if (~isempty(centerConePurityNonOpponentMcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'nonOpponentOnly')))
        scatter(ax, centerConePurityNonOpponentMcenters, surroundConePurityNonOpponentMcenters, markerSize^2, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.8 0.8 0.8]);
    end

    % 2L:1M cone mosaic
    plot(ax, [0 1], 2/3*[1 1], 'k:', 'LineWidth', ff.lineWidth);
    plot(ax, 2/3*[1 1], [0 1], 'k:', 'LineWidth', ff.lineWidth);
    plot(ax, [0 1], [0 1], 'k--', 'LineWidth', ff.lineWidth);

    xlabel(ax, sprintf('center cone purity'));
    ylabel(ax, sprintf('surround cone purity'));

    axis(ax,'equal')
    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1, ...
        'XTickLabel', {'0', '', '.2', '', '.4', '', '.6', '', '.8', '', '1'}, ...
        'YTickLabel', {'0', '', '.2', '', '.4', '', '.6', '', '.8', '', '1'});
    
    title(ax, plotTitle);
    PublicationReadyPlotLib.applyFormat(ax,ff);
    drawnow;


    % The histograms
    centerPurities = [];
    surroundPurities = [];

    if (~isempty(centerConePurityOpponentLcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'opponentOnly')))
        centerPurities = cat(1, centerPurities(:), centerConePurityOpponentLcenters(:));
        surroundPurities = cat(1, surroundPurities(:), surroundConePurityOpponentLcenters(:));
    end
    
    if (~isempty(centerConePurityOpponentMcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'opponentOnly')))
        centerPurities = cat(1, centerPurities(:), centerConePurityOpponentMcenters(:));
        surroundPurities = cat(1, surroundPurities(:), surroundConePurityOpponentMcenters(:));
    end

    if (~isempty(centerConePurityNonOpponentLcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'nonOpponentOnly')))
        centerPurities = cat(1, centerPurities(:), centerConePurityNonOpponentLcenters(:));
        surroundPurities = cat(1, surroundPurities(:), surroundConePurityNonOpponentLcenters(:));
    end

    if (~isempty(centerConePurityNonOpponentMcenters)) && ((strcmp(whichCells, 'all')) || (strcmp(whichCells, 'nonOpponentOnly')))
        centerPurities = cat(1, centerPurities(:), centerConePurityNonOpponentMcenters(:));
        surroundPurities = cat(1, surroundPurities(:), surroundConePurityNonOpponentMcenters(:));
    end

    conePurityBinWidth = 1/25;

    % The surround purity histogram
    hFigSurroundPurityHistogram = figure(figNo+1); clf;
    
    % Change the width to be short
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    ff.figureSize(1) = 550;

    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFigSurroundPurityHistogram,ff);
    ax = theAxes{1,1};

    histogram(ax, surroundPurities, ...
        'FaceColor', [0.4 0.4 0.4], ...
        'FaceAlpha', 0.75, ...
        'EdgeAlpha', 0.0, ...
        'BinWidth', conePurityBinWidth, ...
        'DisplayStyle', 'bar', ...
        'Normalization', 'probability', ...
        'Orientation','horizontal');
    hold(ax, 'on');
    histogram(ax, surroundPurities, ...
        'FaceColor', 'none', ...
        'EdgeColor', 0.15*[1 1 1], ...
        'LineWidth', 1.0, ...
        'BinWidth', conePurityBinWidth, ...
        'DisplayStyle', 'stairs', ...
        'Normalization', 'probability', ...
        'Orientation','horizontal');
    hold(ax, 'off');
    
    ylabel(ax, 'surround cone purity (L/L+M)');
    xlabel(ax, 'probability');

    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'YLim', [0 1],  'YTick', 0:0.1:1, ...
        'YTickLabel', {'0', '', '.2', '', '.4', '', '.6', '', '.8', '', '1'}, ...
        'XLim', [0 1], ...
        'XTick', 0:0.1:1, 'XTickLabel', {'0', '.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9', '1'});
    
    axis(ax, 'square');

    title(ax, plotTitle);
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'XColor', 'none', 'YColor', 'none');
    drawnow;


    % The center purity histogram
    hFigCenterPurityHistogram = figure(figNo+2); clf;

    % Change the height to be short
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    ff.figureSize(1) = 593;
    ff.figureSize(2) = 500;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFigCenterPurityHistogram,ff);
    ax = theAxes{1,1};

    histogram(ax, centerPurities,  ...
        'FaceColor', [0.4 0.4 0.4], ...
        'FaceAlpha', 0.75, ...
        'EdgeAlpha', 0.0, ...
        'BinWidth', conePurityBinWidth, ...
        'DisplayStyle', 'bar', ...
        'Normalization', 'probability', ...
        'Orientation','vertical');
    hold(ax, 'on');
    histogram(ax, centerPurities, ...
        'FaceColor', 'none', ...
        'EdgeColor', 0.15*[1 1 1], ...
        'LineWidth', 1.0, ...
        'BinWidth', conePurityBinWidth, ...
        'DisplayStyle', 'stairs', ...
        'Normalization', 'probability', ...
        'Orientation','vertical');
    hold(ax, 'off');


    xlabel(ax, 'center cone purity (L/L+M)');
    ylabel(ax, 'probability');

    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [0 1],  'XTick', 0:0.1:1,  ...
        'XTickLabel', {'0', '.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9', '1'}, ...
        'YLim', [0 1], ...
        'YTick', 0:0.1:1, 'YTickLabel', {'0', '', '.2', '', '.4', '', '.6', '', '.8', '', '1'});
    axis(ax, 'square');

    title(ax, plotTitle);
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'XColor', 'none', 'YColor', 'none');
    drawnow;

end


function [hFig, opponentLcenterIndices, nonOpponentLcenterIndices, ...
    opponentMcenterIndices, nonOpponentMcenterIndices] = generateResponseRatioPhasePlot(...
    figNo, ...
    responseStrengthRatiosLcenterDominated, ...
    responseStrengthRatiosMcenterDominated, ...
    phaseDifferenceDegsLcenterDominated, ...
    phaseDifferenceDegsMcenterDominated, ...
    maxResponseStrengthRatio, plotTitle)

    xL = responseStrengthRatiosLcenterDominated .* cosd(phaseDifferenceDegsLcenterDominated);
    yL = responseStrengthRatiosLcenterDominated .* sind(phaseDifferenceDegsLcenterDominated);
    xM = responseStrengthRatiosMcenterDominated .* cosd(phaseDifferenceDegsMcenterDominated);
    yM = responseStrengthRatiosMcenterDominated .* sind(phaseDifferenceDegsMcenterDominated);

    idx = find(...
        (max(abs(xL)) > maxResponseStrengthRatio) | ...
        (max(abs(xM)) > maxResponseStrengthRatio) | ...
        (max(abs(yL)) > maxResponseStrengthRatio) | ...
        (max(abs(yM)) > maxResponseStrengthRatio));
    if (~isempty(idx))
        fprintf(2, 'There are %d cells that plot outside the range of ResponseRangePhasePlot\n', numel(idx));
    end


    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
    hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    if ((numel(responseStrengthRatiosLcenterDominated) == 1) && (numel(responseStrengthRatiosMcenterDominated) == 0)) || ...
       ((numel(responseStrengthRatiosLcenterDominated) == 0) && (numel(responseStrengthRatiosMcenterDominated) == 1))
        markerSize = ff.markerSize;
        lineWidth = ff.lineWidth;
    else
        markerSize = ff.markerSize-8;
        lineWidth = ff.lineWidth/2;
    end

    % Render the polar grid
    theRadii = 0.2:0.2:1.0;
    theAngles = 0:30:360;
    renderPolarGrid(ax, theRadii, theAngles, maxResponseStrengthRatio, ff.lineWidth/2, [0.4 0.4 0.4], plotTitle);



    opponentLcenterIndices = find(xL<0);
    nonOpponentLcenterIndices = find(xL>=0);
    opponentMcenterIndices = find(xM<0);
    nonOpponentMcenterIndices = find(xM>=0);

    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    % The opponent L-center mRGCs
    scatter(ax, xL(opponentLcenterIndices)/maxResponseStrengthRatio, yL(opponentLcenterIndices)/maxResponseStrengthRatio, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor',  LconeSTFcolor*0.5, 'MarkerFaceColor',  LconeSTFcolor)

    % The non-opponent L-center mRGCs
    scatter(ax, xL(nonOpponentLcenterIndices)/maxResponseStrengthRatio, yL(nonOpponentLcenterIndices)/maxResponseStrengthRatio, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);
    
    % The opponent M-center mRGCs
    scatter(ax, xM(opponentMcenterIndices)/maxResponseStrengthRatio, yM(opponentMcenterIndices)/maxResponseStrengthRatio, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ...
        'MarkerEdgeColor',  MconeSTFcolor*0.5, 'MarkerFaceColor',  MconeSTFcolor)
    
    % The non-opponent M-center mRGCs
    scatter(ax, xM(nonOpponentMcenterIndices)/maxResponseStrengthRatio, yM(nonOpponentMcenterIndices)/maxResponseStrengthRatio, markerSize^2, ...
        'LineWidth', lineWidth, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7, ....
        'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.75 0.75 0.75]);

    
    axis(ax,'equal')
    set(ax, 'FontSize', ff.axisFontSize);
    set(ax, 'XLim', [-1.02 1.02], 'YLim', [-1.02 1.02], 'XTick', [], 'YTick', []);
    PublicationReadyPlotLib.applyFormat(ax,ff);
    set(ax, 'XColor', 'none', 'YColor', 'none');

    drawnow;
end


function renderPolarGrid(ax, theRadii, theAngles, maxR, lineWidth, lineColor, plotTitle)
    hold(ax, 'on')
    for idx = 1:numel(theRadii)
        radius = theRadii(idx);
        plot(ax,radius*cosd(0:1:360), radius*sind(0:1:360), 'k-', 'Color', lineColor, 'LineWidth', lineWidth);
    end

    for idx = 1:numel(theAngles)
        angle = theAngles(idx);
        plot(ax, [0 cosd(angle)], [0 sind(angle)], 'k-', 'Color', lineColor, 'LineWidth', lineWidth);
    end

    title(ax, sprintf('%s', plotTitle));
end

function renderPolarPlot(ax, theRadii, theAngles, markerSizes, lineWidth, markerFaceAlpha, markerEdgeAlpha, markerFaceColor, markerEdgeColor)
    xL = theRadii .* cosd(theAngles);
    yL = theRadii .* sind(theAngles);
    plot(ax, xL,yL, 'k-');

    for iPoint = 1:numel(xL)
        scatter(ax, xL(iPoint),yL(iPoint), (markerSizes(iPoint))^2, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', markerFaceAlpha, 'MarkerEdgeAlpha', markerEdgeAlpha, ...
            'MarkerEdgeColor', markerFaceColor, 'MarkerFaceColor', markerEdgeColor);
    end

end

function hFig = visualizeSingleCellConeWeightsMapsAndConeIsolatingSTFs(theMRGCMosaic, iRGC, centerConeDominance, ...
    sfSupport, theLconeIsolatingSTF, theMconeIsolatingSTF, unwrapSTFphase, employNativeSTFphaseUnwrapMethod, ...
    forceBandPassSTFphaseUnwrap, unwrapLMopponentSTFphase)
    
    hFig = figure(1000); clf;
	set(hFig, 'Position', [10 10 1200 1000], 'Color', [1 1 1]);
    plotWidth = 0.42;
    plotHeight = 0.40;
	axConeWeightsMap = axes('Position', [0.07 0.09 plotWidth plotHeight]);
	axConeWeightsLineWeightingFunctions = axes('Position', [0.07 0.56 plotWidth plotHeight]);
	axConeIsolatingSTFamplitudes = axes('Position', [0.57 0.56 plotWidth plotHeight]);
    axConeIsolatingSTFphases = axes('Position', [0.57 0.09 plotWidth plotHeight]);

    theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(iRGC,:);
    [scaleBarDegs, scaleBarMicrons, ...
	 spatialSupportTickSeparationArcMin, ...
	 spatialSupportCenterDegs, ...
		 domainVisualizationLimits, ...
		 domainVisualizationTicks, ...
		 domainVisualizationLimitsSingleRF, ...
		 domainVisualizationTicksSingleRF] = ...
		 	RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);

	
    domainVisualizationTicksSingleRF.x = domainVisualizationTicksSingleRF.x(1:2:end);
    figureFormat = PublicationReadyPlotLib.figureComponents('1x1 standard figure'); 
	figNo = [];
    [centerLineWeightingFunctions, surroundLineWeightingFunctions] = ...
     	RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(figNo, ...
                theMRGCMosaic, iRGC, '', ...
                'domainVisualizationLimits', domainVisualizationLimitsSingleRF, ...
                'domainVisualizationTicks', domainVisualizationTicksSingleRF, ...
                'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
                'fixedScaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'noGrid', true, ...
                'plotTitle', '', ...
                'figureHandle', hFig, ...
                'axesHandle', axConeWeightsMap, ...
                'figureFormat', figureFormat, ...
                'renderLineWeightingFunctionPlots', false);
    

    whichMeridian = 'horizontal';
    RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
			'', [], [], ...
			spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
			centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
			'axesToRenderIn', axConeWeightsLineWeightingFunctions, ...
			'domainVisualizationLimits', domainVisualizationLimitsSingleRF(1:2), ...
			'domainVisualizationTicks', domainVisualizationTicksSingleRF.x);

    set(axConeWeightsLineWeightingFunctions, 'XTickLabel', {});
    title(axConeWeightsLineWeightingFunctions, ...
        sprintf('RGC %d/%d at (%2.2f,%2.2f)', iRGC, theMRGCMosaic.rgcsNum, ...
		theRGCpositionDegs(1), theRGCpositionDegs(2)));

    xlabel( axConeWeightsLineWeightingFunctions, '');


    % The STF amplitudes
    LconeSTFcolor = RGCMosaicConstructor.constants.LcenterColor;
    MconeSTFcolor = RGCMosaicConstructor.constants.McenterColor;

    theLconeIsolatingSTFmaxAmplitude = max(theLconeIsolatingSTF.amplitudeSpectrum(:));
    theMconeIsolatingSTFmaxAmplitude = max(theMconeIsolatingSTF.amplitudeSpectrum(:));

    maxAmplitude = max([theLconeIsolatingSTFmaxAmplitude theMconeIsolatingSTFmaxAmplitude]);
    plot(axConeIsolatingSTFamplitudes, sfSupport, theLconeIsolatingSTF.amplitudeSpectrum/maxAmplitude, ...
            '-', 'LineWidth', figureFormat.lineWidth+2, 'Color', [0 0 0]);
    hold(axConeIsolatingSTFamplitudes,'on');
    plot(axConeIsolatingSTFamplitudes, sfSupport, theLconeIsolatingSTF.amplitudeSpectrum/maxAmplitude, ...
            '-', 'LineWidth', figureFormat.lineWidth, 'Color', LconeSTFcolor);

    plot(axConeIsolatingSTFamplitudes, sfSupport, theMconeIsolatingSTF.amplitudeSpectrum/maxAmplitude, ...
            '-', 'Color', MconeSTFcolor, 'LineWidth', figureFormat.lineWidth+2);
    plot(axConeIsolatingSTFamplitudes, sfSupport, theMconeIsolatingSTF.amplitudeSpectrum/maxAmplitude, ...
            '-', 'Color', MconeSTFcolor, 'LineWidth', figureFormat.lineWidth);

    for iSF = 1:numel(sfSupport)
        plot(axConeIsolatingSTFamplitudes, sfSupport(iSF), theLconeIsolatingSTF.amplitudeSpectrum(iSF)/maxAmplitude, ...
            'o', 'LineWidth', figureFormat.lineWidth/2, ...
            'MarkerFaceColor', LconeSTFcolor, ...
            'MarkerEdgeColor', LconeSTFcolor*0.5, ...
            'MarkerSize', figureFormat.markerSize);
    end

    for iSF = 1:numel(sfSupport)
        plot(axConeIsolatingSTFamplitudes, sfSupport(iSF), theMconeIsolatingSTF.amplitudeSpectrum(iSF)/maxAmplitude, ...
            'o',  'LineWidth', figureFormat.lineWidth/2, ...
            'MarkerFaceColor', MconeSTFcolor, ...
            'MarkerEdgeColor', MconeSTFcolor*0.5, ...
            'MarkerSize', figureFormat.markerSize);
    end
  

    axis(axConeIsolatingSTFamplitudes, 'square');
    xLims = [0.008 100];
    xTicks = [0.01 0.1 1 10 100];
    xTickLabels = {'.01', '.1', '1', '10', '100'};
    yLims = [1e-3*5 1.05];
    set(axConeIsolatingSTFamplitudes, ...
        'XScale', 'log', 'XLim', xLims, ...
        'XTick', xTicks, 'XTickLabel', xTickLabels, ...
        'YScale', 'log', 'YLim', yLims, 'YTick', [0.01 0.1 1], 'YTickLabel', {'.01', '.1', '1'});
    if (centerConeDominance == cMosaic.LCONE_ID)
        title(axConeIsolatingSTFamplitudes, 'L-cone center')
    else
        title(axConeIsolatingSTFamplitudes, 'M-cone center')
    end
    xlabel(axConeIsolatingSTFamplitudes, '')
    ylabel(axConeIsolatingSTFamplitudes, 'amplitude')
    PublicationReadyPlotLib.applyFormat(axConeIsolatingSTFamplitudes,figureFormat);
    %PublicationReadyPlotLib.offsetAxes(axConeIsolatingSTFamplitudes, figureFormat, xLims, yLims);


    if (unwrapSTFphase)
        if (employNativeSTFphaseUnwrapMethod)
            jumpTolerangeRadians = 1.2 * pi;
            LconePhaseDegs = unwrap(theLconeIsolatingSTF.phaseSpectrum/180*pi, jumpTolerangeRadians)/pi*180;
            MconePhaseDegs = unwrap(theMconeIsolatingSTF.phaseSpectrum/180*pi, jumpTolerangeRadians)/pi*180;
        else
            theComplexLconeIsolatingSTF = theLconeIsolatingSTF.amplitudeSpectrum .* exp(1i*theLconeIsolatingSTF.phaseSpectrum/180*pi);
            LconePhaseDegs = RGCMosaicAnalyzer.compute.unwrappedPhaseOfBandPassSTF(theComplexLconeIsolatingSTF,'Deg', forceBandPassSTFphaseUnwrap);
        
            theComplexMconeIsolatingSTF = theMconeIsolatingSTF.amplitudeSpectrum .* exp(1i*theMconeIsolatingSTF.phaseSpectrum/180*pi);
            MconePhaseDegs = RGCMosaicAnalyzer.compute.unwrappedPhaseOfBandPassSTF(theComplexMconeIsolatingSTF,'Deg', forceBandPassSTFphaseUnwrap);
        end
    else
        LconePhaseDegs = theLconeIsolatingSTF.phaseSpectrum;
        MconePhaseDegs = theMconeIsolatingSTF.phaseSpectrum;
    end

    % Make sure the lowest SF phase in in the range [-360: 360]
    if (LconePhaseDegs(1) > 360)
        LconePhaseDegs = LconePhaseDegs - 360;
    elseif (LconePhaseDegs(1) < -360)
        LconePhaseDegs = LconePhaseDegs + 360;
    end

    if (MconePhaseDegs(1) > 360)
        MconePhaseDegs = MconePhaseDegs - 360;
    elseif (MconePhaseDegs(1) < -360)
        MconePhaseDegs = MconePhaseDegs + 360;
    end


    if (MconePhaseDegs(1) < -180) && (LconePhaseDegs(1) > MconePhaseDegs(1) + 330)
       MconePhaseDegs = MconePhaseDegs + 360; 
    end

    if (MconePhaseDegs(1) > 180) && (LconePhaseDegs(1) < MconePhaseDegs(1) - 330)
        MconePhaseDegs = MconePhaseDegs - 360; 
    end

    if (LconePhaseDegs(1) < -180) && (MconePhaseDegs(1) > LconePhaseDegs(1) + 330)
       LconePhaseDegs = LconePhaseDegs + 360; 
    end

    if (LconePhaseDegs(1) > 180) && (MconePhaseDegs(1) < LconePhaseDegs(1) - 330)
        LconePhaseDegs = LconePhaseDegs - 360; 
    end





    plot(axConeIsolatingSTFphases, sfSupport, LconePhaseDegs, '-', ...
        'LineWidth', figureFormat.lineWidth+2, 'Color', [0 0 0]);
    hold(axConeIsolatingSTFphases, 'on');
    plot(axConeIsolatingSTFphases, sfSupport, LconePhaseDegs, '-', ...
        'LineWidth', figureFormat.lineWidth, 'Color', LconeSTFcolor);

    plot(axConeIsolatingSTFphases, sfSupport, MconePhaseDegs, '-', ...
        'LineWidth', figureFormat.lineWidth+2, 'Color', [0 0 0]);
    plot(axConeIsolatingSTFphases, sfSupport, MconePhaseDegs, '-', ...
        'LineWidth', figureFormat.lineWidth, 'Color', MconeSTFcolor);

    for iSF = 1:numel(sfSupport)
       plot(axConeIsolatingSTFphases, sfSupport(iSF), LconePhaseDegs(iSF), ...
                'ko-', 'LineWidth', figureFormat.lineWidth/2, ...
                'MarkerFaceColor', LconeSTFcolor, 'MarkerEdgeColor', [0 0 0], ...
                'MarkerSize', figureFormat.markerSize);
    end
    
    for iSF = 1:numel(sfSupport)
       plot(axConeIsolatingSTFphases, sfSupport(iSF), MconePhaseDegs(iSF), ...
                'ko-', 'LineWidth', figureFormat.lineWidth/2, ...
                'MarkerFaceColor', MconeSTFcolor, 'MarkerEdgeColor', [0 0 0],...
                'MarkerSize', figureFormat.markerSize);
    end
    
    deltaPhaseDegs = LconePhaseDegs-MconePhaseDegs;

    if (unwrapLMopponentSTFphase)
        if (employNativeSTFphaseUnwrapMethod)
            jumpTolerangeRadians = 1.2 * pi;
            deltaPhaseDegs = unwrap(deltaPhaseDegs/180*pi, jumpTolerangeRadians)/pi*180;
        else
            theComplexDeltaSTF = ...
                theLconeIsolatingSTF.amplitudeSpectrum .* exp(1i*theLconeIsolatingSTF.phaseSpectrum/180*pi) - ...
                theMconeIsolatingSTF.amplitudeSpectrum .* exp(1i*theMconeIsolatingSTF.phaseSpectrum/180*pi);

            deltaPhaseDegs = RGCMosaicAnalyzer.compute.unwrappedPhaseOfBandPassSTF(theComplexDeltaSTF,'Deg', forceBandPassSTFphaseUnwrap);
        end 
    end


    plot(axConeIsolatingSTFphases, sfSupport, deltaPhaseDegs, ...
                'k--', 'LineWidth', figureFormat.lineWidth)
    
    axis(axConeIsolatingSTFphases, 'square');
    yLims = [-360-180 360+180];
    set(axConeIsolatingSTFphases, ...
        'XScale', 'log', 'XLim', xLims, ...
        'XTick', xTicks, 'XTickLabel', xTickLabels, ...
        'YScale', 'linear', 'YLim', yLims, 'YTick', -720:90:720);
    xlabel(axConeIsolatingSTFphases, 'spatial frequency (c/deg)')
    ylabel(axConeIsolatingSTFphases, 'L-M phase (degs)');

    drawnow;
    

	PublicationReadyPlotLib.applyFormat(axConeIsolatingSTFphases,figureFormat);
    %PublicationReadyPlotLib.offsetAxes(axConeIsolatingSTFphases, figureFormat, xLims, yLims);
    
end



function [theLconeIsolatingSTF, theMconeIsolatingSTF, sfSupport, centerConeDominance] = extractTheAnalyzedSTFslices(theMRGCMosaic, theRGCindex, stimParams, ...
            theMRGCMosaicResponseTemporalSupportSeconds, ...
            theLconeIsolatingSpatioTemporalResponses2DSTF, ...
			theMconeIsolatingSpatioTemporalResponses2DSTF, ...
            fixedOptimalOrientation, ...
            incrementExcitatoryLowSFresponsePhaseDegs)

    
    % Method to return cone connectivity stats for theSubregion ('center' or 'surround') oftheRGCindex
    s = theMRGCMosaic.singleCellConnectivityStats(theRGCindex, 'center', ...
        'minConeWeightIncluded', 0.01, ...
        'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);

    visualizeFullAndMaximalExcursionSTF = false;
    
    theLconeIsolatingSTF = struct();
    theMconeIsolatingSTF = struct();
    sfSupport = stimParams.spatialFrequencyCPD;
    centerConeDominance = s.dominantConeType;

    if (centerConeDominance == cMosaic.LCONE_ID)

        % L-center dominated cell.
        % Get the L-cone STF slice at the requested fixedOptimalOrientation
        % and save the actual orientation so we can use it below to obtain
        % the M-cone STF
        [~, theCenterConeIsolatingOptimalOrientation, ...
	     theLconeIsolatingSTF.phaseSpectrum, ...
	     theLconeIsolatingSTF.amplitudeSpectrum, ...
	     theLconeIsolatingSTF.fullAmplitudeSpectra, ...
	     theLconeIsolatingSTF.fullPhaseSpectra] = RGCMosaicConstructor.helper.simulateExperiment.maximalExcursionSTFfrom2DSTF(...
				    stimParams.orientationDegs, stimParams.spatialFrequencyCPD, stimParams.spatialPhasesDegs, ...
				    theMRGCMosaicResponseTemporalSupportSeconds, ...
				    theLconeIsolatingSpatioTemporalResponses2DSTF, ...
				    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
				    'fixedOptimalOrientation', fixedOptimalOrientation, ...
				    'axFullSTF', [], ...
    		        'axSTFslice', []);
         
        % Get the M-cone STF slice at the same orientation as the L-cone STF
        [~, ~, ...
	      theMconeIsolatingSTF.phaseSpectrum, ...
	      theMconeIsolatingSTF.amplitudeSpectrum, ...
	      theMconeIsolatingSTF.fullAmplitudeSpectra, ...
	      theMconeIsolatingSTF.fullPhaseSpectra] = RGCMosaicConstructor.helper.simulateExperiment.maximalExcursionSTFfrom2DSTF(...
				    stimParams.orientationDegs, stimParams.spatialFrequencyCPD, stimParams.spatialPhasesDegs, ...
				    theMRGCMosaicResponseTemporalSupportSeconds, ...
				    theMconeIsolatingSpatioTemporalResponses2DSTF, ...
				    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
				    'fixedOptimalOrientation', theCenterConeIsolatingOptimalOrientation, ...
				    'axFullSTF', [], ...
    		        'axSTFslice', []);
    else
        % M-center dominated cell.
        % Get the M-cone STF slice at the requested fixedOptimalOrientation
        % and save the actual orientation so we can use it below to obtain
        % the L-cone STF
        [~, theCenterConeIsolatingOptimalOrientation, ...
	     theMconeIsolatingSTF.phaseSpectrum, ...
	     theMconeIsolatingSTF.amplitudeSpectrum, ...
	     theMconeIsolatingSTF.fullAmplitudeSpectra, ...
	     theMconeIsolatingSTF.fullPhaseSpectra] = RGCMosaicConstructor.helper.simulateExperiment.maximalExcursionSTFfrom2DSTF(...
				    stimParams.orientationDegs, stimParams.spatialFrequencyCPD, stimParams.spatialPhasesDegs, ...
				    theMRGCMosaicResponseTemporalSupportSeconds, ...
				    theMconeIsolatingSpatioTemporalResponses2DSTF, ...
				    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
				    'fixedOptimalOrientation', fixedOptimalOrientation, ...
				    'axFullSTF', [], ...
    		        'axSTFslice', []);
         
        % Get the L-cone STF slice at the same orientation as the M-cone STF
        [~, ~, ...
	      theLconeIsolatingSTF.phaseSpectrum, ...
	      theLconeIsolatingSTF.amplitudeSpectrum, ...
	      theLconeIsolatingSTF.fullAmplitudeSpectra, ...
	      theLconeIsolatingSTF.fullPhaseSpectra] = RGCMosaicConstructor.helper.simulateExperiment.maximalExcursionSTFfrom2DSTF(...
				    stimParams.orientationDegs, stimParams.spatialFrequencyCPD, stimParams.spatialPhasesDegs, ...
				    theMRGCMosaicResponseTemporalSupportSeconds, ...
				    theLconeIsolatingSpatioTemporalResponses2DSTF, ...
				    'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
				    'fixedOptimalOrientation', theCenterConeIsolatingOptimalOrientation, ...
				    'axFullSTF', [], ...
    		        'axSTFslice', []);
    end

    theLconeIsolatingSTF.phaseSpectrum = theLconeIsolatingSTF.phaseSpectrum - incrementExcitatoryLowSFresponsePhaseDegs;
    theLconeIsolatingSTF.fullPhaseSpectra = theLconeIsolatingSTF.fullPhaseSpectra - incrementExcitatoryLowSFresponsePhaseDegs;
    theMconeIsolatingSTF.phaseSpectrum = theMconeIsolatingSTF.phaseSpectrum - incrementExcitatoryLowSFresponsePhaseDegs;
    theMconeIsolatingSTF.fullPhaseSpectra = theMconeIsolatingSTF.fullPhaseSpectra - incrementExcitatoryLowSFresponsePhaseDegs;
end
