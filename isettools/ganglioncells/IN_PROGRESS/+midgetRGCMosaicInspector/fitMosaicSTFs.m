function fitMosaicSTFs(mosaicCenterParams, rfModelParams, opticsParams,  maxRGCsNum)

    midgetRGCMosaicInspector.say('Fitting mosaic STFs using the difference of Gaussians model');

    % Generate the frozen mosaic filename
    frozenMosaicFileName = midgetRGCMosaicInspector.frozenMosaicFileName(...
        mosaicCenterParams, rfModelParams.H1cellIndex, opticsParams);
    
    % Load the frozen midget RGC mosaic
    load(frozenMosaicFileName, 'theMidgetRGCmosaic');

    % Ask the user to specify the optics position for which responses were saved
    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        opticsPositionDegs = input('\nEnter the optics position that was used to compute the responses ([x y]): ');
    end

    % Generate the responses filename
    responsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        frozenMosaicFileName, opticsPositionDegs);
       

    % Load the mosaic responses
    load(responsesFileName, 'theMidgetRGCMosaicResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');


    % Examine every 15 degs
    deltaAngle = 22.5;
    theMeridianAngles  = 0:deltaAngle:(180-deltaAngle);
    theMeridianRadius = 1.5 * sqrt(2.0);

    % Compute RGC indices to fit along each meridian
    RGCindicesToFit = cell(1, numel(theMeridianAngles));
    for iMeridianAngle = 1:numel(theMeridianAngles)
        theMeridianAngle = theMeridianAngles(iMeridianAngle);
        idx = midgetRGCMosaicInspector.rgcIndicesAlongMeridianWithAngle(...
            theMeridianAngle, theMeridianRadius, ...
            theMidgetRGCmosaic.rgcRFpositionsDegs, maxRGCsNum);
        RGCindicesToFit{iMeridianAngle} = idx;
    end

    
    visualizeFits = true;
    if (visualizeFits)
        hFigRcDegs = []; hFigRsRcRatios = []; hFigSCintSensRatios = [];
    end


    % Fit STFs along all meridians
    theMeridianFits = cell(1, numel(theMeridianAngles));
    for iMeridianAngle = 1:numel(theMeridianAngles)

       fprintf('Fitting RGCs along the %d deg meridian\n', theMeridianAngles(iMeridianAngle));
       rgcIndicesAlongThisMeridian = RGCindicesToFit{iMeridianAngle};

       fitsAlongThisMeridian = struct();
       fitsAlongThisMeridian.rgcIndicesAlongThisMeridian = rgcIndicesAlongThisMeridian;

       % Fit all RGCs located at this meridian
       [fitsAlongThisMeridian.fittedParams, fitsAlongThisMeridian.fittedSTFs] = fitSelectSTFs(...
            rgcIndicesAlongThisMeridian, ...
            spatialFrequenciesTested, orientationsTested, ...
            theMidgetRGCmosaic, theMidgetRGCMosaicResponses);
    
       
       if (visualizeFits)
            [hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios] = ...
                midgetRGCMosaicInspector.renderSTFfitPlots(...
                    hFigRcDegs, hFigRsRcRatios, hFigSCintSensRatios, ...
                    theMidgetRGCmosaic, ...
                    theMeridianAngles, iMeridianAngle, ...
                    fitsAlongThisMeridian);
       end

        theMeridianFits{iMeridianAngle} = fitsAlongThisMeridian;
    end


    % Append the fittedSTFs structs.
    % Here we use  matfile which lets you read and write to part of variables in a mat file.
    % In that way, if 'theMeridianAngles', 'theMeridianFits' do not exist
    % it will write them to the file, whereas if they do exist it will
    % replace them with the new values
    m = matfile(responsesFileName, 'Writable', true);
    m.theMeridianAngles = theMeridianAngles;
    m.theMeridianFits = theMeridianFits;
    fprintf('Appended theMeridianFits to the responses file: ''%s''.', responsesFileName);
    clear 'm';
end


function [d, fittedSTFs] = fitSelectSTFs(rgcIndicesToAnalyze, ...
    spatialFrequenciesTested, orientationsTested, ...
    theMidgetRGCmosaic, theMidgetRGCMosaicResponses)

    % Allocate memory
    temporalEquivalentEccDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRcDegs = zeros(1, numel(rgcIndicesToAnalyze));
    achievedRsToRcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedKsToKcRatios = zeros(1, numel(rgcIndicesToAnalyze));
    achievedSCintSensRatios = zeros(1, numel(rgcIndicesToAnalyze)); 
    conesNumPooledByTheRFcenter = zeros(1, numel(rgcIndicesToAnalyze)); 
    fittedSTFs = cell(1, numel(rgcIndicesToAnalyze));

    parfor iRGC = 1:numel(rgcIndicesToAnalyze)
        % Target RGC
        theRGCindex = rgcIndicesToAnalyze(iRGC);
        %fprintf('Fitting RGC %d of %d\n', iRGC, numel(rgcIndicesToAnalyze));

        % Do the fitting
        [fittedSTFs{iRGC}, ...
         temporalEquivalentEccDegs(iRGC), ...
         conesNumPooledByTheRFcenter(iRGC), ...
         achievedRcDegs(iRGC), ...
         achievedRsToRcRatios(iRGC), ...
         achievedKsToKcRatios(iRGC), ...
         achievedSCintSensRatios(iRGC)] = fitSingleRGC(...
                theMidgetRGCmosaic, ...
                theMidgetRGCMosaicResponses, ...
                theRGCindex, ...
                spatialFrequenciesTested, ...
                orientationsTested);

    end % iRGC

    d = struct(...
        'temporalEquivalentEccDegs', temporalEquivalentEccDegs, ...
        'conesNumPooledByTheRFcenter', conesNumPooledByTheRFcenter, ...
        'achievedRcDegs', achievedRcDegs, ...
        'achievedRsToRcRatios', achievedRsToRcRatios, ...
        'achievedKsToKcRatios', achievedKsToKcRatios, ...
        'achievedSCintSensRatios', achievedSCintSensRatios ...
        );
end


function [fittedSTF, temporalEquivalentEccDegs, conesNumPooledByTheRFcenter, ...
          achievedRcDegs, achievedRsToRcRatio, achievedKsToKcRatio, achievedSCintSensRatio] = ...
            fitSingleRGC(theMidgetRGCmosaic, theMidgetRGCMosaicResponses, theRGCindex, ...
            spatialFrequenciesTested, orientationsTested)
   
        temporalEquivEccDegs = theMidgetRGCmosaic.temporalEquivalentEccentricityForEccentricity(theMidgetRGCmosaic.rgcRFpositionsDegs(theRGCindex,:));
        temporalEquivalentEccDegs = sqrt(sum(temporalEquivEccDegs.^2,2));

        % Obtain the responses
        theSingleMidgetRGCResponses = squeeze(theMidgetRGCMosaicResponses(:, :, :, theRGCindex));

        % Select STF to fit (orientation that results in max extension into high SFs)
        [theMeasuredSTF, allMeasuredSTFs] = highestExtensionSTF(theSingleMidgetRGCResponses, spatialFrequenciesTested, orientationsTested);

        % Fit the DoG model to the measured STF
        multiStartsNum = 128;
        [RcDegsEstimate, conesNumPooledByTheRFcenter] = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex);
         
        
        [theFittedDoGmodelParams, theDoGmodelFitOfTheMeasuredSTF] = ...
            RTVF.fitDoGmodelToMeasuredSTF(spatialFrequenciesTested, ...
                    theMeasuredSTF, ...
                    RcDegsEstimate, ...
                    [], ...
                    multiStartsNum);


        % The RcDegs
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RcDegs'));
        achievedRcDegs = theFittedDoGmodelParams.finalValues(idx);

        % The Rs/Rc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'RsToRc'));
        achievedRsToRcRatio = theFittedDoGmodelParams.finalValues(idx);

        % The Ks/Kc ratio
        idx = find(strcmp(theFittedDoGmodelParams.names, 'kS/kC'));
        achievedKsToKcRatio = theFittedDoGmodelParams.finalValues(idx);

        % The S/C int sens ratio
        achievedSCintSensRatio = achievedKsToKcRatio * (achievedRsToRcRatio)^2;
        

        fittedSTF = struct();
        fittedSTF.targetRGC = theRGCindex;
        fittedSTF.spatialFrequencySupport = spatialFrequenciesTested;
        fittedSTF.orientationsTested = orientationsTested;
        fittedSTF.allMeasuredSTFs =  allMeasuredSTFs;
        fittedSTF.theMeasuredSTF = theMeasuredSTF;
        fittedSTF.theDoGmodelFitOfTheMeasuredSTF = theDoGmodelFitOfTheMeasuredSTF;
        fittedSTF.theFittedDoGmodelParams = theFittedDoGmodelParams;  
end

function [RcDegs, conesNumPooledByTheRFcenter] = retinalRFcenterRcDegsInitialEstimate(theMidgetRGCmosaic, theRGCindex)

    connectivityVector = full(squeeze(theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix(:, theRGCindex)));
    indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);

    conesNumPooledByTheRFcenter = numel(indicesOfCenterCones);
    coneRcDegs = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones)) * ...
                 theMidgetRGCmosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
    RcDegs = sqrt(conesNumPooledByTheRFcenter)*coneRcDegs;
end

function [theHighestExtensionSTF, theMeasuredSTFs] = highestExtensionSTF(theMidgetRGCMosaicResponses, spatialFrequenciesTested, orientationsTested)

    % Allocate memory
    theMeasuredSTFs = zeros(numel(orientationsTested),numel(spatialFrequenciesTested));

    for iSF = 1:numel(spatialFrequenciesTested)
        for iOri = 1:numel(orientationsTested)
            % Retrieve the mRGC response time-series
            theResponseTimeSeries = squeeze(theMidgetRGCMosaicResponses(iOri, iSF, :));
    
            % Compute the response modulation for this SF
            theMeasuredSTFs(iOri, iSF) = max(theResponseTimeSeries)-min(theResponseTimeSeries);
        end % iORI
    end % iSF

    theMeasuredSTFs = theMeasuredSTFs / max(theMeasuredSTFs(:));

    % Determine the orientation that maximizes the STF extension to high spatial frequencies
    maxSF = nan(1,numel(orientationsTested));
    spatialFrequenciesInterpolated = linspace(spatialFrequenciesTested(1),spatialFrequenciesTested(end), 50);

    for iOri = 1:numel(orientationsTested)
        % Find spatial frequency at which STF drops to 20% of max
        theSTFatThisOri = squeeze(theMeasuredSTFs(iOri,:));
        theSTFatThisOriInterpolated = interp1(spatialFrequenciesTested, theSTFatThisOri, spatialFrequenciesInterpolated);
        [mag, iSFpeak] = max(theSTFatThisOri);
        thresholdSTF = mag * 0.2;

        ii = iSFpeak;
        keepGoing = true; iStop = [];
        while (ii < numel(spatialFrequenciesInterpolated)-1)&&(keepGoing)
            ii = ii + 1;
            if (theSTFatThisOriInterpolated(ii)>=thresholdSTF) && (theSTFatThisOriInterpolated(ii+1)<thresholdSTF)
                keepGoing = false;
                iStop = ii;
            end
        end % while
        if (~isempty(iStop))
            maxSF(iOri) = spatialFrequenciesInterpolated(iStop);
        end
    end % iOri

    % Best orientation
    if (any(isnan(maxSF)))
        theSTFatTheHighestSF = squeeze(theMeasuredSTFs(:,end));
        [~, iBestOri] = max(theSTFatTheHighestSF(:));
    else
        [~, iBestOri] = max(maxSF);
    end
    theHighestExtensionSTF = squeeze(theMeasuredSTFs(iBestOri,:));
end