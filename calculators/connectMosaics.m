function connectMosaics()
   
    % Select mosaics to load
    whichEye = 'right';
    mosaicFOVDegs = 15; %30;
    eccentricitySamplesNumCones = 32; %48;
    eccentricitySamplesNumRGC = 32;
    
    % Connect mosaics only within a central region to save compute time
    connectivityRadiusDeg = 2.5;
    
    % Load data for the analyzed region
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios] = ...
        loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, eccentricitySamplesNumRGC, connectivityRadiusDeg);
    
    % Compute connection matrix between the 2 mosaics
    connectionMatrix = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, conesToRGCratios);

    % Visualized connectivity in a chosen region of interest
    roi.centerDeg = [0 0];
    roi.sizeDeg = [0.1 0.1];
    figNo = 1;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, roi)

    roi.centerDeg = [1 0];
    roi.sizeDeg = [0.25 0.25];
    figNo = 2;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, roi)

    roi.centerDeg = [2 0];
    roi.sizeDeg = [0.5 0.5];
    figNo = 3;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, roi)


end


function visualizeConnectivity(figNo, connectionMatrix, RGCRFPositions, RGCRFSpacings, conePositions, coneSpacings, roi)

    % Convert to microns
    micronsPerDeg = 300;
    roi.center = roi.centerDeg * micronsPerDeg;
    roi.size = roi.sizeDeg * micronsPerDeg;
    

    d = 2*abs(bsxfun(@minus, RGCRFPositions, roi.center));
    mRGCindices = find((d(:,1) <= roi.size(1)) & (d(:,2) <= roi.size(2)));
    
    d = 2*abs(bsxfun(@minus, conePositions, roi.center));
    coneIndices = find((d(:,1) <= roi.size(1)) & (d(:,2) <= roi.size(2)));
    
    % Sampling for contours
    deltaX = 0.1;
    xAxis = (roi.center(1)-roi.size(1)/2): deltaX: (roi.center(1)+roi.size(1)/2);
    yAxis = (roi.center(2)-roi.size(2)/2): deltaX: (roi.center(2)+roi.size(2)/2);
    [X,Y] = meshgrid(xAxis,yAxis);
    
    figure(figNo); clf;
    hold on;
    
    % plot all cones
    plot(conePositions(coneIndices,1), conePositions(coneIndices,2), 'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
    
    % Plot mRGCs
    for m = 1:numel(mRGCindices)
        mRGCindex = mRGCindices(m);
        
        % Find which cones are connected to this RGC
        cIndices = find(connectionMatrix(mRGCindex,:)>0);
 
        % Generate RGC outline based on its cone inputs
        coneDiam = 0.7*mean(coneSpacings(cIndices));
        coneSigma = coneDiam/6;
        
        zLevels = exp([-3 -2 -1]);
        [xRGCEnsembleOutline, yRGCEnsembleOutline, xPeak, yPeak] = generateRGCRFoutlineBasedOnItsConnectivity(...
            connectionMatrix(mRGCindex,cIndices), conePositions(cIndices,:), coneSigma, X, Y, xAxis, yAxis,  zLevels);
        
        % Plot original RGC RF 
        %mPos = RGCRFPositions(mRGCindex,:);
        %xRGCoutline = mPos(1)+cosd(0:20:360)*RGCRFSpacings(mRGCindex)/4;
        %yRGCoutline = mPos(2)+sind(0:20:360)*RGCRFSpacings(mRGCindex)/4;
        %plot(xRGCoutline, yRGCoutline, 'b--', 'LineWidth', 1); 
        
        % PLot RGC RF based on its cone inputs
        whichLevelsToContour = [1];
        multipleContourPlot(xRGCEnsembleOutline, yRGCEnsembleOutline, whichLevelsToContour);
        
        
        for c = 1:numel(cIndices)
            cIndex = cIndices(c);
            connectionStrength = connectionMatrix(mRGCindex,cIndex);
            if (connectionStrength>0)
                plot([conePositions(cIndex,1) xPeak], [conePositions(cIndex,2) yPeak], 'r.-', 'LineWidth', 1); %connectionStrength*3);
                
                if (connectionStrength == 0.0123456789)
                    plot(conePositions(cIndex,1), conePositions(cIndex,2), 'r.', 'MarkerSize', 10);
    
                    fprintf(2,'Cone %d connected to mRGC %d with %f strength\n', cIndex, mRGCindex, connectionStrength);
                    pause
                else
                    fprintf('Cone %d connected to mRGC %d with %f strength\n', cIndex, mRGCindex, connectionStrength);
                end
            end
        end
        
        if (numel(cIndices)>1)
            fprintf(2,'More than 2 inputs\n');
        end
        
        xTickIncrDeg = 0.05;
        xTickIncrMicrons = xTickIncrDeg * micronsPerDeg;
        xTicks = [xTickIncrMicrons:xTickIncrMicrons:1000];
        xTicks = [-fliplr(xTicks) 0 xTicks];
        xTickDegLabel = sprintf('%2.2f\n', xTicks/micronsPerDeg);
        
        axis 'equal'
        lessSpaceMicrons = 4;
        xLims = roi.center(1) + (roi.size(1)/2-lessSpaceMicrons)*[-1 1];
        yLims = roi.center(2) + (roi.size(2)/2-lessSpaceMicrons)*[-1 1];
        set(gca, 'FontSize', 12, 'XLim', xLims, 'YLim', yLims,  ...
            'xTick', xTicks, 'YTick', xTicks, ...
            'xTickLabel', xTickDegLabel, 'yTickLabel', xTickDegLabel);
        xlabel('\it space (deg)');
        box on
        
        drawnow
    end
end

function  [xRGCEnsembleOutline, yRGCEnsembleOutline, xPeak, yPeak] = ...
            generateRGCRFoutlineBasedOnItsConnectivity(...
            connectionsVector, conePositions, coneSigma, X,Y, xAxis, yAxis,  zLevels)
        

    ensembleRF = [];  
    for inputConeIndex = 1:size(conePositions)
        cP = squeeze(conePositions(inputConeIndex,:));
        % raise gain to <1 to better visualize secondary inputs
        gain = double(connectionsVector(inputConeIndex));
        % Make flat top coneProfile
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile = coneProfile.^0.35;
        coneProfile = coneProfile/max(coneProfile(:));
        
        coneRF =  gain * coneProfile;
        if (isempty(ensembleRF))
            ensembleRF = coneRF;
        else
            ensembleRF = ensembleRF + coneRF;
        end
    end
    
    [maxRF, idx] = max(ensembleRF(:));
    [row,col] = ind2sub(size(ensembleRF), idx);
    xPeak = xAxis(col);
    yPeak = yAxis(row);
    
    ensembleRF = ensembleRF / maxRF;
    C = contourc(xAxis, yAxis,ensembleRF, zLevels);
    k = 1;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        for kLevel = 1:numel(zLevels)
            if (level == zLevels(kLevel))
                xRGCEnsembleOutline.level{kLevel} = C(1,k+(1:points));
                yRGCEnsembleOutline.level{kLevel} = C(2,k+(1:points));
            end
        end
        k = k+points+1;
    end
        
end


function multipleContourPlot(xRGCEnsembleOutline, yRGCEnsembleOutline, whichLevels)
    
    for iLevel = 1:numel(whichLevels)
        theLevel = whichLevels(iLevel);
        f = 1:numel(xRGCEnsembleOutline.level{theLevel});
        v = [xRGCEnsembleOutline.level{theLevel}(:) yRGCEnsembleOutline.level{theLevel}(:)];
        patch('Faces', f, 'Vertices', v, 'FaceColor', [0.8 0.8 0.8], ...
            'FaceAlpha', 0.5, 'EdgeColor', [0 0 0.6], 'LineWidth', 1.0);
    end
    
end

function [RGCRFPositions, RGCRFSpacings, conePositions, coneSpacings, conesToRGCratios] = ...
        loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, eccentricitySamplesNumRGC, analyzedRadiusDeg)
    
    micronsPerDeg = 300;
    workingRadiusMicrons = analyzedRadiusDeg*micronsPerDeg;
    
    % Load the positions of the mRGC mosaic
    neuronalType = 'mRGC';
    RGCRFPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumRGC);
    
    % Compute spacing and cone-to-RGC ratios for each  RGCRF position
    [RGCRFSpacings, conesToRGCratios] = mRGCStats(RGCRFPositions, 128, whichEye);

    % Load the positions of the cone mosaic
    neuronalType = 'cone';
    conePositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumCones);
    
    % Compute cone spacing for each cone location
    coneSpacings = coneStats(conePositions, 128, whichEye);
    
    % Only keep cones  within the working radius
    eccCones = sqrt(sum(conePositions.^2,2));
    idx = find(eccCones<workingRadiusMicrons);
    conePositions = conePositions(idx,:);
    coneSpacings = coneSpacings(idx);
    
    % Only keep RGC within the working radius
    eccRGC = sqrt(sum(RGCRFPositions.^2,2));
    idx = find(eccRGC<workingRadiusMicrons);
    RGCRFPositions = RGCRFPositions(idx,:);
    RGCRFSpacings = RGCRFSpacings(idx);
    conesToRGCratios = conesToRGCratios(idx);
    
    fprintf('Loaded %2.0f RGCs\n', size(RGCRFPositions,1));
    fprintf('Loaded %2.0f cones\n', size(conePositions,1));
end


function connectionMatrix = computeConnectionMatrix(RGCRFPositions, conePositions, RGCRFSpacings, conesToRGCratios)

    conesNum = size(conePositions,1);
    RGCsNum = size(RGCRFPositions,1);
    
    connectionMatrix = zeros(RGCsNum, conesNum, 'single');
    
    RGCinputsNum = zeros(RGCsNum,1);
    
    coneEcc = sqrt(sum(conePositions.^2,2));
    [~,sortedConeIndices] = sort(coneEcc, 'ascend');
    
        
    % Connect every cone to its nearest RGC 
    for cIndex = 1:conesNum
        theTargetConeIndex = sortedConeIndices(cIndex);
        theTargetConePos = conePositions(theTargetConeIndex,:);
        
        % distance to all RGCs
        distanceToAllRGCs = sqrt(sum((bsxfun(@minus, RGCRFPositions, theTargetConePos)).^2,2));
        
        % find the closest RGC 
        [~,closestRGCindex] = min(distanceToAllRGCs);
        
        % Choose the spacing of that RGC as the search radius
        searchRadius = RGCRFSpacings(closestRGCindex);
        
        % Find nearby RGCs (within a searchRadius from the target cone)
        possibleRGCindicesForConnection = find(distanceToAllRGCs <= searchRadius);
          
        % Examine which of those nearby RGCs to connect to
        if (~isempty(possibleRGCindicesForConnection))

            % Likelihood of connecting is inversely proportional to distance
            willAcceptAnotherConnection = 1./distanceToAllRGCs(possibleRGCindicesForConnection);
            willAcceptAnotherConnection(willAcceptAnotherConnection>100) = 100;
            
            % Do not connect connections to RGCs that already have a connection 
            for k = 1:numel(possibleRGCindicesForConnection)
                theRGCindex = possibleRGCindicesForConnection(k);
                if (RGCinputsNum(theRGCindex) > 0)
                    willAcceptAnotherConnection(k) = 0;
                end
            end
            
            % Find RGC with P of connecting to
            [maxAcceptance,idx] = max(willAcceptAnotherConnection);
            
            if (maxAcceptance > 0)
                matchedRGCindex = possibleRGCindicesForConnection(idx);
                %fprintf('Connecting 1st cone to RGC %d\n',  matchedRGCindex);
                connectionStrength = 1;
            else
                % All possible RGCs have at least 1 input, so find the RGC
                % with the minimal cone connections that is also closest to the target cone
                % The min number of cone inputs in this collection of RGCs
                minConeInputs = min(RGCinputsNum(possibleRGCindicesForConnection));
                
                % The RGCs which have this minimum number of inputs
                idx = find(RGCinputsNum(possibleRGCindicesForConnection) == minConeInputs);
                RGCindicesWithLeastConeInputs = possibleRGCindicesForConnection(idx);
                
                % Sort them according to their distance to the target cone
                [~,idx] = min(distanceToAllRGCs(RGCindicesWithLeastConeInputs));
                
                % Select the RGC with the closest distance
                matchedRGCindex = RGCindicesWithLeastConeInputs(idx);
                
                % Comute connection weight
                meanExpectedInputs = conesToRGCratios(matchedRGCindex);
                currentInputs = minConeInputs;
                powerF = 2.0;
                connectionStrength = exp(-(currentInputs/meanExpectedInputs)^powerF);
                
                if (minConeInputs == 1)
                    postfix = 'nd';
                elseif (minConeInputs == 2)
                    postfix = 'rd';
                else
                    postfix = 'th';
                end

                if (minConeInputs >= 2)
                    fprintf(2, 'Connecting %d%s cone to RGC %d (which has mean cone inputs of %2.3f) with strength:%f\n', ...
                        RGCinputsNum(matchedRGCindex)+1,postfix,  matchedRGCindex, conesToRGCratios(matchedRGCindex), connectionStrength);
                else
                    %fprintf('Connecting %d%s cone to RGC %d (which has mean cone inputs of %2.3f) with strength:%f\n', ...
                    %    RGCinputsNum(matchedRGCindex)+1,postfix,  matchedRGCindex, conesToRGCratios(matchedRGCindex), connectionStrength);
                end
                
                
            end
            
            connectionMatrix(matchedRGCindex,theTargetConeIndex) = single(connectionStrength);
                 
            % Change this RGC position to match that of the cone
            %RGCRFPositionsNew(matchedRGCindex,:) = theTargetConePos;

            % Update the input connections to this RGC
            RGCinputsNum(matchedRGCindex) = RGCinputsNum(matchedRGCindex)+1;
        end
        
    end % cIndex

    %
    netConeWeights = sum(connectionMatrix,2);
    connectionMatrix = bsxfun(@times, connectionMatrix, 1./netConeWeights);
end




function rfPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum)
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    
    
    RGCMosaicFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_prctile%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, 99));

    load(RGCMosaicFileName, 'rfPositionsHistory');
    rfPositions = double(squeeze(rfPositionsHistory(end,:,:)));   
end


function coneSpacingInMicrons = coneStats(rfPositions, nSamples, whichEye)
    % Find range of retinal positions in microns that we need to compute density for
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    idx = find(eccentricitiesInMicrons<2.1);
    t = rfPositions(idx,:);
    tt = pdist2(t,t);
    tt(tt==0) = Inf;
    minSeparationMicrons = min(tt(:));
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Support of retinal positions in mm
    xPosMM = [0 logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples)]/1e3;
    
    switch whichEye
        case 'left'
            theView = 'left eye retina';
        case 'right'
            theView = 'right eye retina';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end

    WatsonRGCModelObj = WatsonRGCModel();
    
    % Convert retinal mm to visual degs
    eccDegs = WatsonRGCModelObj.rhoMMsToDegs(xPosMM);
    
    % Compute cone density map
    [coneDensity2DMap, coneMeridianDensities, densitySupportMM, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DConeRFDensity(eccDegs, theView);
    
     % Make sure results are returned in retinal mm units
    assert((strcmp(densityUnits, WatsonRGCModelObj.retinalMMDensityUnits)) && ...
           (strcmp(supportUnits, WatsonRGCModelObj.retinalMMEccUnits)), ...
           sprintf('Expected mm units, but got ''%s'' and ''%s'' instead.', supportUnits, densityUnits)); 
       
    % Convert the density map into a spacing map
    coneSpacing2DMapMM = WatsonRGCModelObj.spacingFromDensity(coneDensity2DMap);
    
    % Convert to microns from mm
    coneSpacing2DMapMicrons = coneSpacing2DMapMM*1e3;
    densitySupportMicrons = densitySupportMM*1e3;

    % Create a scatterred interpolant function  so we can
    % compute mRGC RF spacing at the actual mRGC positions
    [X,Y] = meshgrid(squeeze(densitySupportMicrons(1,:)), squeeze(densitySupportMicrons(2,:)));
    Fspacing = scatteredInterpolant(X(:),Y(:),coneSpacing2DMapMicrons(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    coneSpacingInMicrons = Fspacing(rfPositions(:,1), rfPositions(:,2));
    
end

function [mRGCSpacingInMicrons, conesToRGCratio, eccentricitiesInMicrons] = mRGCStats(rfPositions, nSamples, whichEye)

    % Find range of retinal positions in microns that we need to compute density for
    eccentricitiesInMicrons = sqrt(sum(rfPositions .^ 2, 2));
    idx = find(eccentricitiesInMicrons<2.1);
    t = rfPositions(idx,:);
    tt = pdist2(t,t);
    tt(tt==0) = Inf;
    minSeparationMicrons = min(tt(:));
    maxEccMicrons = max(eccentricitiesInMicrons(:));
    
    % Support of retinal positions in mm
    xPosMM = [0 logspace(log10(minSeparationMicrons), log10(maxEccMicrons+minSeparationMicrons), nSamples)]/1e3;
    
    switch whichEye
        case 'left'
            theView = 'left eye retina';
        case 'right'
            theView = 'right eye retina';
        otherwise
            error('Which eye must be either ''left'' or ''right'', not ''%s''.', whichEye)
    end

    WatsonRGCModelObj = WatsonRGCModel();
    
    % Convert retinal mm to visual degs
    eccDegs = WatsonRGCModelObj.rhoMMsToDegs(xPosMM);

    % Compute mRGC density map
    [mRGCDensity2DMap, mRGCMeridianDensities, densitySupportMM, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DmRGCRFDensity(eccDegs, theView);

    % Density is for both types of mRGCs (ON + OFF), so we need density for
    % one type, which is half (assuming equal numerosities of ON and OFF cells)
    mRGCDensity2DMap = 0.5*mRGCDensity2DMap;
    
    % Compute cone density map
    [coneDensity2DMap, coneMeridianDensities, densitySupport, ...
        horizontalMeridianLabel, verticalMeridianLabel, densityLabel, ...
        supportUnits, densityUnits] = WatsonRGCModelObj.compute2DConeRFDensity(eccDegs, theView);
   
    
    % Make sure results are returned in retinal mm units
    assert((strcmp(densityUnits, WatsonRGCModelObj.retinalMMDensityUnits)) && ...
           (strcmp(supportUnits, WatsonRGCModelObj.retinalMMEccUnits)), ...
           sprintf('Expected mm units, but got ''%s'' and ''%s'' instead.', supportUnits, densityUnits)); 

    % Compute cone to mRGC ratio map
    conesToMRGCratio2Dmap = coneDensity2DMap./mRGCDensity2DMap;
    conesToMRGCratio2Dmap(conesToMRGCratio2Dmap<1) = 1;
    
    % Convert the density map into a spacing map
    mRGCSpacing2DMapMM = WatsonRGCModelObj.spacingFromDensity(mRGCDensity2DMap);
    
    % Convert to microns from mm
    mRGCSpacing2DMapMicrons = mRGCSpacing2DMapMM*1e3;
    densitySupportMicrons = densitySupportMM*1e3;

    % Create a scatterred interpolant function  so we can
    % compute mRGC RF spacing at the actual mRGC positions
    [X,Y] = meshgrid(squeeze(densitySupportMicrons(1,:)), squeeze(densitySupportMicrons(2,:)));
    Fspacing = scatteredInterpolant(X(:),Y(:),mRGCSpacing2DMapMicrons(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    mRGCSpacingInMicrons = Fspacing(rfPositions(:,1), rfPositions(:,2));
    
    % Create a scatterred interpolant function  so we can
    % compute cone-to-mRGC ratios at the actual mRGC positions
    Fratio = scatteredInterpolant(X(:),Y(:),conesToMRGCratio2Dmap(:), 'linear');

    % Evaluate the interpolant function at the requested rfPositions
    conesToRGCratio = Fratio(rfPositions(:,1), rfPositions(:,2));
    
end



