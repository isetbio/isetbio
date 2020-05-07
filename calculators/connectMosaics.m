function connectMosaics()
   
    % Select mosaics to load
    whichEye = 'right';
    mosaicFOVDegs = 15;
    eccentricitySamplesNumCones = 32;  
    eccentricitySamplesNumRGC = 32; 
    maxMovementPercentileCones = 20;
    maxMovementPercentileRGC = 20;
    bestIterationCones = Inf;
    bestIterationRGC = 95;
    
    
    % Connect mosaics only within a central region to save compute time
    connectivityRadiusDeg = 2.5;
    
    % Load data for the analyzed region
    [RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios] = ...
        loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, eccentricitySamplesNumRGC, ...
        maxMovementPercentileCones, maxMovementPercentileRGC, ...
         bestIterationCones,  bestIterationRGC, connectivityRadiusDeg);
    
    % Compute connection matrix between the 2 mosaics
    connectionMatrix = computeConnectionMatrix(RGCRFPositionsMicrons, conePositionsMicrons, RGCRFSpacingsMicrons, conesToRGCratios);

    % Visualized connectivity in a chosen region of interest
    roi.centerDeg = [0 0];
    roi.sizeDeg = [0.07 0.07];
    figNo = 1;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios, roi)

    roi.centerDeg = [1.0 0];
    roi.sizeDeg = [0.15 0.15];
    figNo = 2;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios, roi)

    
    roi.centerDeg = [2.0 0];
    roi.sizeDeg = [0.25 0.25];
    figNo = 2;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios, roi)

    roi.centerDeg = [4.7 0];
    roi.sizeDeg = [0.3 0.3];
    figNo = 3;
    visualizeConnectivity(figNo, connectionMatrix, RGCRFPositionsMicrons, RGCRFSpacingsMicrons, conePositionsMicrons, coneSpacingsMicrons, conesToRGCratios, roi)

    

end


function visualizeConnectivity(figNo, connectionMatrix, RGCRFPositions, RGCRFSpacings, conePositions, coneSpacings, conesToRGCratios, roi)

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
    subplot(1,2,1);
    plot(conePositions(coneIndices,1), conePositions(coneIndices,2), 'ro', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 10); hold on;
    plot(RGCRFPositions(mRGCindices,1), RGCRFPositions(mRGCindices,2), 'ko',  'MarkerSize', 16, 'LineWidth', 1.0);
    legend({'cones', 'mRGC'});
    drawnow;
%     for m = 1:numel(mRGCindices)
%         mPos = RGCRFPositions(mRGCindices(m),:);
%         connectionStrengths = full(connectionMatrix(mRGCindices(m),coneIndices));
%         idx = find(connectionStrengths>0);
%         for c = 1:numel(idx)
%             cPos =  conePositions(coneIndices(idx),:);
%             plot([mPos(1) cPos(1)], [mPos(2) cPos(2)], 'r-');
%             drawnow;
%         end
%     end
    coneToRGCRatio = numel(coneIndices)/numel(mRGCindices)
    title(sprintf('actual cone to RGC ratio: %2.3f', coneToRGCRatio));
    
                
    subplot(1,2,2);
    hold on;
    
    % plot all cones
    plotCones = false;
    if (plotCones)
        plot(conePositions(coneIndices,1), conePositions(coneIndices,2), 'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
    end
    
    % Plot mRGCs
    for m = 1:numel(mRGCindices)
        mRGCindex = mRGCindices(m);
        
        % Find which cones are connected to this RGC
        cIndices = find(connectionMatrix(mRGCindex,:)>0);
 
        % Generate RGC outline based on its cone inputs
        coneDiam = 0.7*mean(coneSpacings(cIndices));
        coneSigma = coneDiam/6;
        
        zLevels = exp([-3 -1]);
        [xRGCEnsembleOutline, yRGCEnsembleOutline, xPeak, yPeak] = generateRGCRFoutlineBasedOnItsConnectivity(...
            full(connectionMatrix(mRGCindex,cIndices)), conePositions(cIndices,:), coneSigma, X, Y, xAxis, yAxis, zLevels);
        
        % Plot original RGC RF 
        %mPos = RGCRFPositions(mRGCindex,:);
        %xRGCoutline = mPos(1)+cosd(0:20:360)*RGCRFSpacings(mRGCindex)/4;
        %yRGCoutline = mPos(2)+sind(0:20:360)*RGCRFSpacings(mRGCindex)/4;
        %plot(xRGCoutline, yRGCoutline, 'b--', 'LineWidth', 1); 
        
        % Plot RGC RF based on its cone inputs
        whichLevelsToContour = [1 2];
        multipleContourPlot(xRGCEnsembleOutline, yRGCEnsembleOutline, whichLevelsToContour, xPeak, yPeak, conesToRGCratios(mRGCindex));
        
        
        for c = 1:numel(cIndices)
            cIndex = cIndices(c);
            connectionStrength = full(connectionMatrix(mRGCindex,cIndex));
            if (connectionStrength>0)
                plot([conePositions(cIndex,1) xPeak], [conePositions(cIndex,2) yPeak], 'r-', 'LineWidth', 1); %connectionStrength*3);
                plot(conePositions(cIndex,1), conePositions(cIndex,2), 'ro', 'MarkerFaceColor', [1 0.4 0.4], 'LineWidth', 1);
                
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
        lessSpaceMicrons = 0;
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
        coneProfile = coneProfile.^0.3;
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


function multipleContourPlot(xRGCEnsembleOutline, yRGCEnsembleOutline, whichLevels, xPeak, yPeak, conesToRGCratio)
    
    for iLevel = 1:numel(whichLevels)
        theLevel = whichLevels(iLevel);
        maxLevel = max(whichLevels);
        f = 1:numel(xRGCEnsembleOutline.level{theLevel});
        v = [xRGCEnsembleOutline.level{theLevel}(:) yRGCEnsembleOutline.level{theLevel}(:)];
        patch('Faces', f, 'Vertices', v, 'FaceColor', (1-0.8*theLevel/maxLevel)*[0.8 0.8 0.8], ...
            'FaceAlpha', 0.5, 'EdgeColor', [0 0 0.6], 'EdgeAlpha', 0, 'LineWidth', 1.0);
    end
    text(xPeak-2, yPeak-1, sprintf('%2.2f', conesToRGCratio), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
end

function [RGCRFPositions, RGCRFSpacings, conePositions, coneSpacings, conesToRGCratios] = ...
        loadData(whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, eccentricitySamplesNumRGC, ...
        maxMovementPercentileCones, maxMovementPercentileRGC,  bestIterationCones,  bestIterationRGC, analyzedRadiusDeg)
    
    micronsPerDeg = 300;
    workingRadiusMicrons = analyzedRadiusDeg*micronsPerDeg;
    
    % Load the positions of the mRGC mosaic
    neuronalType = 'mRGC';
    RGCRFPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumRGC, maxMovementPercentileRGC,  bestIterationRGC);
    
    % Compute spacing and cone-to-RGC ratios for each  RGCRF position
    [RGCRFSpacings, conesToRGCratios] = mRGCStats(RGCRFPositions, 128, whichEye);

    % Load the positions of the cone mosaic
    neuronalType = 'cone';
    conePositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNumCones, maxMovementPercentileCones,  bestIterationCones);
    
    % Compute cone spacing for each cone location
    coneSpacings = coneStats(conePositions, 128, whichEye);
    
    % Only keep cones  within the working radius
    eccCones = sqrt(sum(conePositions.^2,2));
    idx = find((eccCones<workingRadiusMicrons) & (conePositions(:,1)>=-20));
    conePositions = conePositions(idx,:);
    coneSpacings = coneSpacings(idx);
    
    % Only keep RGC within the working radius
    eccRGC = sqrt(sum(RGCRFPositions.^2,2));
    idx = find((eccRGC<workingRadiusMicrons) & (RGCRFPositions(:,1)>=-20));
    RGCRFPositions = RGCRFPositions(idx,:);
    RGCRFSpacings = RGCRFSpacings(idx);
    conesToRGCratios = conesToRGCratios(idx);
    
    minRGCEccDegs = min(RGCRFPositions(:,1))/micronsPerDeg;
    maxRGCEccDegs = max(RGCRFPositions(:,1))/micronsPerDeg;
    
    minRGCEccCones = min(conePositions(:,1))/micronsPerDeg;
    maxRGCEccCones = max(conePositions(:,1))/micronsPerDeg;
    
    fprintf('Loaded %2.0f RGCs from %2.2f-%2.2f degs\n', size(RGCRFPositions,1), minRGCEccDegs, maxRGCEccDegs);
    fprintf('Loaded %2.0f cones from %2.2f-%2.2f degs\n', size(conePositions,1), minRGCEccCones, maxRGCEccCones);
end


function connectionMatrix = computeConnectionMatrix(RGCRFPositions, conePositions, RGCRFSpacings, conesToRGCratios)

    conesNum = size(conePositions,1);
    RGCsNum = size(RGCRFPositions,1);
    
    fprintf('Connecting %d cones to %d RGCs\n', conesNum, RGCsNum);
    pause
    connectionMatrix = sparse(RGCsNum, conesNum) ; % zeros(RGCsNum, conesNum, 'single');
    
    RGCinputsNum = zeros(RGCsNum,1);
    
    % Sort cones according to their eccentricity so we can start connecting
    % from fovea towards the periphery
    coneEcc = sqrt(sum(conePositions.^2,2));
    [~,sortedConeIndices] = sort(coneEcc, 'ascend');
    
    % Connect every cone to its nearest RGC 
    for cIndex = 1:conesNum
        theTargetConeIndex = sortedConeIndices(cIndex);
        theTargetConePos = conePositions(theTargetConeIndex,:);
        
        % compute distance to all RGCs
        distanceToAllRGCs = sqrt(sum((bsxfun(@minus, RGCRFPositions, theTargetConePos)).^2,2));
        
        % find the closest RGC 
        [~,closestRGCindex] = min(distanceToAllRGCs);
        
        % Set the search radius proportional to the spacing of the closest RGC
        searchRadius = 0.75*RGCRFSpacings(closestRGCindex);
        
        % Find nearby RGCs (within a searchRadius from the target cone)
        possibleRGCindicesForConnection = find(distanceToAllRGCs <= searchRadius);
          
        % Examine which of those nearby RGCs to connect to
        if (~isempty(possibleRGCindicesForConnection))

            % Likelihood of connecting is inversely proportional to distance
            willAcceptAnotherConnection = 1./distanceToAllRGCs(possibleRGCindicesForConnection);
            willAcceptAnotherConnection(willAcceptAnotherConnection>100) = 100;
            
            % Do not connect to RGCs that already have a connection 
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

                if (minConeInputs >= 1)
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

    % Normalize cone weights by dividing with the net cone weight
    netConeWeights = sum(connectionMatrix,2);
    
    % Or Normalize cone weights by dividing with the max cone weight
    maxConeWeights = max(connectionMatrix,[],2);
    
    normalizingFactor = maxConeWeights;
    connectionMatrix = bsxfun(@times, connectionMatrix, 1./normalizingFactor);
end




function rfPositions = getPositions(neuronalType, whichEye, mosaicFOVDegs, eccentricitySamplesNum, maxMovementPercentile, bestIteration)
    % Save filename
    p = getpref('IBIOColorDetect');
    mosaicDir = strrep(p.validationRootDir, 'validations', 'sideprojects/MosaicGenerator'); 
    
    mosaicFileName = fullfile(mosaicDir, sprintf('progress_%s_%s_Mosaic%2.1fdegs_samplesNum%d_maxMovPrctile%d.mat', ...
        whichEye, neuronalType, mosaicFOVDegs, eccentricitySamplesNum, maxMovementPercentile));

    load(mosaicFileName, 'rfPositionsHistory', 'reTriangulationIterations');
    if (~isinf(bestIteration))
        [~, targetIterationIndex] = min(abs(reTriangulationIterations - bestIteration));
    else
        targetIterationIndex = numel(reTriangulationIterations);
    end
    
    fprintf('Loading %s mosaic with %d neurons from iteration %d\n', neuronalType, size(rfPositionsHistory,2), reTriangulationIterations (targetIterationIndex));
    rfPositions = double(squeeze(rfPositionsHistory(targetIterationIndex,:,:)));   
    
    contrastMeridianDensitiesFromActualMosaicsToModel(neuronalType,rfPositions, mosaicFOVDegs);
end

function contrastMeridianDensitiesFromActualMosaicsToModel(neuronalType, rfPositions, mosaicFOVDegs)
    
    meridianDensities = [];

    if (strcmp(neuronalType, 'cone'))
        figNo = 99;
    else
        figNo = 100;
    end
    figure(figNo); clf;
    
    samplesNum = 30;
    
    eccAxis = logspace(log10(0.01), log10(mosaicFOVDegs/2), samplesNum);
    
    subplot(1,2,1);
    targetPositionsDegs(:,1) = eccAxis;
    targetPositionsDegs(:,2) = eccAxis*0;
    [targetHorizRFPositionsMicrons, meanLocalSpacingMicrons] = computeMeanLocalSpacing(rfPositions, targetPositionsDegs);
    plot(WatsonRGCModel.rhoMMsToDegs(1e-3*targetHorizRFPositionsMicrons(:,1)), WatsonRGCModel.densityFromSpacing(meanLocalSpacingMicrons*1e-3), 'rs-', 'LineWidth', 1.0); hold on
    
    targetPositionsDegs(:,1) = eccAxis * 0;
    targetPositionsDegs(:,2) = eccAxis;
    [targetHorizRFPositionsMicrons, meanLocalSpacingMicrons] = computeMeanLocalSpacing(rfPositions, targetPositionsDegs);
    plot(WatsonRGCModel.rhoMMsToDegs(1e-3*targetHorizRFPositionsMicrons(:,2)), WatsonRGCModel.densityFromSpacing(meanLocalSpacingMicrons*1e-3), 'bs-', 'LineWidth', 1.0);
    
    targetPositionsDegs(:,1) = -fliplr(eccAxis);
    targetPositionsDegs(:,2) = eccAxis * 0;
    [targetHorizRFPositionsMicrons, meanLocalSpacingMicrons] = computeMeanLocalSpacing(rfPositions, targetPositionsDegs);
    plot(-WatsonRGCModel.rhoMMsToDegs(1e-3*targetHorizRFPositionsMicrons(:,1)), WatsonRGCModel.densityFromSpacing(meanLocalSpacingMicrons*1e-3), 'rs-', 'LineWidth', 1.0);
    
    targetPositionsDegs(:,1) = eccAxis * 0;
    targetPositionsDegs(:,2) = -fliplr(eccAxis);
    [targetHorizRFPositionsMicrons, meanLocalSpacingMicrons] = computeMeanLocalSpacing(rfPositions, targetPositionsDegs);
    plot(-WatsonRGCModel.rhoMMsToDegs(1e-3*targetHorizRFPositionsMicrons(:,2)), WatsonRGCModel.densityFromSpacing(meanLocalSpacingMicrons*1e-3), 'bs-', 'LineWidth', 1.0);
    xlabel('ecc(degs)');
    ylabel(sprintf('density from mosaic (%s per mm^2)', neuronalType));
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.1 10], 'YLim', [1e3 1e5*3], 'XTick', [0.1 0.3 1 3 10], 'YTick', [1000 3000 10000 30000  100000]);
    grid on;
    
    obj = WatsonRGCModel();
    subplot(1,2,2);
    rightEyeVisualFieldMeridianName = obj.enumeratedMeridianNames{1};
    eccUnits = obj.visualDegsEccUnits;
    densityUnits = obj.retinalMMDensityUnits;
    if (strcmp(neuronalType, 'cone'))
        [coneRFSpacing, density, rightEyeRetinalMeridianName] = obj.coneRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
    else
        [mRGCRFSpacing, density, rightEyeRetinalMeridianName] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
         % half density because the reported density is for both ON and OFF mRGCs
         density = density /2 ;
    end
    
    plot(eccAxis, density, 'rs-', 'LineWidth', 1.0); hold on
    rightEyeVisualFieldMeridianName = obj.enumeratedMeridianNames{2};
    if (strcmp(neuronalType, 'cone'))
        [coneRFSpacing, density, rightEyeRetinalMeridianName] = obj.coneRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
    else
        [mRGCRFSpacing, density, rightEyeRetinalMeridianName] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
         % half density because the reported density is for both ON and OFF mRGCs
         density = density /2;
    end
    plot(eccAxis, density, 'gs-', 'LineWidth', 1.0);
    
    rightEyeVisualFieldMeridianName = obj.enumeratedMeridianNames{3};
    if (strcmp(neuronalType, 'cone'))
        [coneRFSpacing, density, rightEyeRetinalMeridianName] = obj.coneRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
    else
        [mRGCRFSpacing, density, rightEyeRetinalMeridianName] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
         % half density because the reported density is for both ON and OFF mRGCs
         density = density /2 ;
    end
    plot(eccAxis, density, 'bs-', 'LineWidth', 1.0);
    
    rightEyeVisualFieldMeridianName = obj.enumeratedMeridianNames{4};
    if (strcmp(neuronalType, 'cone'))
        [coneRFSpacing, density, rightEyeRetinalMeridianName] = obj.coneRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
    else
        [mRGCRFSpacing, density, rightEyeRetinalMeridianName] = obj.mRGCRFSpacingAndDensityAlongMeridian(eccAxis, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        % half density because the reported density is for both ON and OFF mRGCs
        density = density /2 ;
    end
    plot(eccAxis, density, 'ms-', 'LineWidth', 1.0);
    xlabel('ecc(degs)');
    ylabel(sprintf('density from model (%s per mm^2)', neuronalType));
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.1 10], 'YLim', [1e3 1e5*3], 'XTick', [0.1 0.3 1 3 10], 'YTick', [1000 3000 10000 30000 100000]);
    grid on
    drawnow
end

function [targetRFPositionsMicrons, meanLocalSpacingMicrons] = computeMeanLocalSpacing(rfPositions, targetPositionsDegs)
    
    targetPositionsMicrons = 1e3*WatsonRGCModel.rhoMMsToDegs(targetPositionsDegs);
    
    plotNeighbors = false;
    targetsNum = size(targetPositionsDegs,1);
    
    for targetIndex = 1:targetsNum
        d = sum((bsxfun(@minus, rfPositions, targetPositionsMicrons(targetIndex,:))).^2,2);
        [~,idx(targetIndex)] = min(d);
    end
    targetRFPositionsMicrons = rfPositions(idx,:);
    [D,I] = pdist2(rfPositions, targetRFPositionsMicrons, 'euclidean', 'Smallest', 7);
    
    meanLocalSpacingMicrons = zeros(targetsNum,1);
    if (plotNeighbors)
        figure(100); clf;
    end
    
    for targetIndex = 1:targetsNum
        for neigborIndex = 2:size(I,1)
            rfIndex = I(neigborIndex,targetIndex);
            if (plotNeighbors)
                plot(rfPositions(rfIndex,1), rfPositions(rfIndex,2), 'ko',  'MarkerSize', 8);
                hold on;
            end
        end
        meanLocalSpacingMicrons(targetIndex) = mean(D(:,targetIndex));
        if (plotNeighbors)
            plot(targetRFPositionsMicrons(targetIndex,1), targetRFPositionsMicrons(targetIndex,2), 'rs', 'MarkerSize', 12);
            plot(targetRFPositionsMicrons(targetIndex,1)+meanLocalSpacingMicrons(targetIndex)*cosd(0:20:360), ...
                targetRFPositionsMicrons(targetIndex,2)+meanLocalSpacingMicrons(targetIndex)*sind(0:20:360), 'r-', 'MarkerSize', 12);
        end
    end
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



