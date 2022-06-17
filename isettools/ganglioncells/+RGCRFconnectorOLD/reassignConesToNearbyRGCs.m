function [RGCRFinputs, RGCRFweights] = reassignConesToNearbyRGCs(...
            operationMode, ...
            theInputConeMosaic, RGCRFinputs, RGCRFweights, ...
            localConeToRGCDensityRatio, ...
            wiringParams, visualizationParams)
% Re-assign cones from one RGC to a nearby RGC with less inputs
%
% Syntax:
%   [RGCRFinputs, RGCRFweights, ] = ...
%        RGCRFconnector.reassignConesToNearbyRGCsWithLessInputs(...
%                    theInputConeMosaic, RGCRFinputs, RGCRFweights, ...
%                    wiringParams, visualizationParams)
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    % Unpack inputs
    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;
    allConeSpacings = theInputConeMosaic.coneRFspacingsMicrons;
    allConeTypes = theInputConeMosaic.coneTypes;

    rgcsNum = numel(RGCRFinputs);

    % Find centroids and # of cone inputs to each RGC.
    [RGCRFcentroids, coneInputsNum] = RGCRFconnector.centroidsFromConeInputs(...
        RGCRFinputs, RGCRFweights, allConePositions);

    % Find spacings of RGCs
    RGCRFspacings = RGCmodels.Watson.convert.positionsToSpacings(RGCRFcentroids);

    groupCounts = zeros(1,100);
    for iRGC = 1:rgcsNum
        coneInputsNum(iRGC) = numel(RGCRFinputs{iRGC});
        theCount = coneInputsNum(iRGC);
        groupCounts(theCount) = groupCounts(theCount)+1;
    end
    fprintf('Max cones/RGC: %d\n', max(coneInputsNum));
    fprintf('Min cones/RGC: %d\n', min(coneInputsNum));


    if ((visualizationParams.exportStagesAsPDF) || (visualizationParams.generateProgressionVideo))
        % PDFfilename
        if (wiringParams.sequentialConeTypeWiring)
            PDFfilename = sprintf('WiringSteps_W%1.2f_SequentialConeTypeWiring_Stage2', wiringParams.chromaticSpatialVarianceTradeoff);
        else
            PDFfilename = sprintf('WiringSteps_W%1.2f_NonSequentialConeTypeWiring_Stage2', wiringParams.chromaticSpatialVarianceTradeoff);
        end

        % Video setup
        visualizationParams.generateProgressionVideo = true;
        if (visualizationParams.generateProgressionVideo)
            % Video setup
            videoOBJ = VideoWriter(PDFfilename, 'MPEG-4');
            videoOBJ.FrameRate = 10;
            videoOBJ.Quality = 100;
            videoOBJ.open();
        else
            videoOBJ = '';
        end

        hFig = figure(100); clf;
        set(hFig, 'Name', 'Stage 2', 'Color', [1 1 1]);
        ax = subplot('Position', [0.1 0.1 0.8 0.8]);
     end

    maxPassesNum = 3;
    iPass = 0;
    successfulTransfersNum = 1;
    
    while (successfulTransfersNum>0)&&(iPass < maxPassesNum)
        iPass = iPass + 1;
        fprintf(2, 'PASS %d\n', iPass);
        if (visualizationParams.generateProgressionVideo)
            stageName = sprintf('Stage-2 (initial)');
            [hFig, ax] = setupFigure(stageName);
            [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
                RGCRFinputs, [], ...
                theInputConeMosaic, ...
                'titleString', stageName, ...
                'superimposeConeInputWiring', true, ...
                'displayRGCID', false, ...
                'pauseAfterEachRGCisRendered', ~true, ...
                'figureHandle', hFig, ...
                'videoOBJ', videoOBJ, ...
                'axesHandle', ax);
        end
        

        % Continue with the most populated RFs and go all the way down to RF with 2 cone inputs.
        for targetConeInputsNum = numel(groupCounts):-1:2
            if (groupCounts(targetConeInputsNum) > 0)
    
                % Find the RGCs that contain this many targetConeInputsNum
                targetRGCindices = find(coneInputsNum == targetConeInputsNum);
                fprintf('\n\nWorking on RGC group with %d cone inputs which contains %d RGCs\n',...
                    targetConeInputsNum, numel(targetRGCindices));
            
                successfulTransfersNum = 0;
    
                for kk = 1:numel(targetRGCindices)
                    targetSourceRGCindex = targetRGCindices(kk);
                    
                    % Find the maxNearbyRGCsNum RGCindices
                    [~, nearbyRGCindices] = pdist2(RGCRFcentroids, RGCRFcentroids(targetSourceRGCindex,:), '', ...
                        'smallest', wiringParams.maxNearbyRGCsNum);
                    
                    % Exclude the targetSourceRGCindex, which has 0 distance
                    nearbyRGCindices = setdiff(nearbyRGCindices, targetSourceRGCindex);
    
                    % Exclude nearbyRGCs that are further than a maxDistance
                    maxDistance = 1.25*RGCRFspacings(targetSourceRGCindex);
                    dd = pdist2(RGCRFcentroids(nearbyRGCindices,:), RGCRFcentroids(targetSourceRGCindex,:));
                    nearbyRGCindices = nearbyRGCindices(find(dd<=maxDistance));
    
                    % Try to reassign a cone input
                    [status, RGCRFinputs, RGCRFweights, RGCRFcentroids, RGCRFspacings, coneInputsNum] = ...
                        RGCRFconnector.tryToReassignConeInput( ...
                            targetSourceRGCindex, ...
                            nearbyRGCindices, ...
                            RGCRFinputs, RGCRFweights, ...
                            RGCRFcentroids, RGCRFspacings, coneInputsNum, ...
                            localConeToRGCDensityRatio, ...
                            wiringParams, allConePositions, allConeSpacings, allConeTypes);
    
                    if (strcmp(status, 'success'))
                        successfulTransfersNum = successfulTransfersNum + 1;
                    end                  
    
                    % Visualize what hapenned   
                    if (visualizationParams.generateProgressionVideo) && (strcmp(status, 'success'))
                        stageName = sprintf('Stage-2/Pass-%d (%d-cone input group, transfer #%d)', iPass, targetConeInputsNum, successfulTransfersNum);
                        [hFig, ax] = setupFigure(stageName);
                        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
                            RGCRFinputs, [], ...
                            theInputConeMosaic, ...
                            'titleString', stageName, ...
                            'superimposeConeInputWiring', true, ...
                            'displayRGCID', false, ...
                            'pauseAfterEachRGCisRendered', ~true, ...
                            'figureHandle', hFig, ...
                            'videoOBJ', videoOBJ, ...
                            'axesHandle', ax);
                     end
    
                end % for kk
    
                fprintf('There were %d successful cone transfers\n', successfulTransfersNum);
                
            end % if
        end % for targetConeInputsNum
    end % while

    iPass = iPass + 1;
    fprintf(2, 'PASS %d. Break RGCs with large inter-input distances. \n', iPass);
    
    successfulTransfersNum = 0;

    for iRGC = 1:rgcsNum
        coneInputs = RGCRFinputs{iRGC};
        conePos = allConePositions(coneInputs,:);
        maxInputConeDistance = RGCRFconnector.maximalInterInputDistance(conePos);

        if (maxInputConeDistance > 1.5*RGCRFspacings(iRGC))

            % Find which cone to transfer: the one most distant
            % from the centroid
            d = pdist2(conePos, RGCRFcentroids(iRGC,:));
            [~, idx] = max(d(:));
            coneIndexToTransfer = coneInputs(idx);

            % Find optimal nearby destination RGC and do the transfer
            [RGCRFinputs, RGCRFweights, RGCRFcentroids, RGCRFspacings] = ...
                RGCRFconnector.transferSpecificConeInputFromSourceRGCToNearbyRGC(...
                RGCRFinputs, RGCRFweights, RGCRFcentroids, ...
                iRGC, coneIndexToTransfer, ...
                localConeToRGCDensityRatio(iRGC), ...
                wiringParams, ...
                allConePositions, allConeSpacings, allConeTypes ...
                );

            % Update transfers num
            successfulTransfersNum = successfulTransfersNum + 1;

            % Visualize what hapenned
            if (visualizationParams.generateProgressionVideo)
                stageName = sprintf('Stage-2/Pass-%d (large inter-cone distance group, transfer #%d)', iPass, successfulTransfersNum);
                [hFig, ax] = setupFigure(stageName);
                [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
                    RGCRFinputs, [], ...
                    theInputConeMosaic, ...
                    'titleString', stageName, ...
                    'superimposeConeInputWiring', true, ...
                    'displayRGCID', false, ...
                    'pauseAfterEachRGCisRendered', ~true, ...
                    'figureHandle', hFig, ...
                    'videoOBJ', videoOBJ, ...
                    'axesHandle', ax);
            end

        end
    end


    if (visualizationParams.exportStagesAsPDF) 
        stageName = 'Stage-2';
        [hFig, ax] = setupFigure(stageName);
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            RGCRFinputs, [], ...
            theInputConeMosaic, ...
            'titleString', stageName, ...
            'superimposeConeInputWiring', true, ...
            'displayRGCID', false, ...
            'pauseAfterEachRGCisRendered', ~true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);

        NicePlot.exportFigToPDF(sprintf('%s.pdf', PDFfilename), hFig, 300);
    
        if (visualizationParams.generateProgressionVideo)
            % Close video stream
            videoOBJ.close();
        end

    end
end

function [hFig, ax] = setupFigure(FigureTitle)
    hFig = figure(100); clf;
    set(hFig, 'Name', FigureTitle, 'Color', [1 1 1], 'Position', [10 10 800 800]);
    ax = subplot('Position', [0.1 0.1 0.8 0.8]);
end