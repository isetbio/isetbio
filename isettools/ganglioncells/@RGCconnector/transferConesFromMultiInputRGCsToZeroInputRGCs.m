function transferConesFromMultiInputRGCsToZeroInputRGCs(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.addParameter('optimizationCenter', 'patchCenter', @(x)(ismember(x, {'patchCenter', 'visualFieldCenter'})));

    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;
    optimizationCenter = p.Results.optimizationCenter;

    % Compute the # of cone inputs for all RGCs
    allRGCsConeInputsNum = squeeze(sum(obj.coneConnectivityMatrix,1));

    % Find out how many RGCs have zero-input
    zeroInputRGCindices = find(allRGCsConeInputsNum == 0);
    fprintf('There are %d RGCs with zero inputs\n', numel(zeroInputRGCindices));

    if (isempty(zeroInputRGCindices))
        return;
    end

    % Video setup
    if (generateProgressVideo)
        videoOBJ = VideoWriter('Step3', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    else
        videoOBJ = [];
    end
    
    % Do this as a function of eccentricity for different groups of RGCs
    coneInputNumerosityGroups = unique(allRGCsConeInputsNum);
    coneInputNumerosityGroups = sort(coneInputNumerosityGroups, 'descend');

    
    for iGroup = 1:numel(coneInputNumerosityGroups)
        if (~isempty(zeroInputRGCindices))
            rgcIndicesForGroup = find(allRGCsConeInputsNum == coneInputNumerosityGroups(iGroup));
    
            targetedRGCCentroids = obj.RGCRFcentroidsFromInputs(rgcIndicesForGroup,:);
            % Compute ecc from centroids
            switch (optimizationCenter)
                case 'visualFieldCenter'
                    ecc = sum(targetedRGCCentroids.^2,2);
                case 'patchCenter'
                    diff = bsxfun(@minus, targetedRGCCentroids, obj.inputConeMosaic.eccentricityMicrons);
                    ecc = sum(diff.^2,2);
            end
        
            % Sort RGCs according to their ecc
            [~, idx] = sort(ecc, 'ascend');
            targetedRGCindices = rgcIndicesForGroup(idx);

            % Do it
            zeroInputRGCindices = doIt(obj, targetedRGCindices, zeroInputRGCindices, allRGCsConeInputsNum, videoOBJ);
        end % (~isempty(zeroInputRGCindices))


    end % iGroup
   

    if (generateProgressVideo)
        videoOBJ.close();
    end
end

function zeroInputRGCindices = doIt(obj, targetedRGCindices, zeroInputRGCindices, allRGCsConeInputsNum, videoOBJ)

    for iRGC = 1:numel(targetedRGCindices)

        if (~isempty(zeroInputRGCindices))
            theSourceRGCindex = targetedRGCindices(iRGC);
            
            % Retrieve the cone inputs
            theSourceRGCinputConeIndices = find(squeeze(obj.coneConnectivityMatrix(:,theSourceRGCindex))>0);
            theSourceRGCinputConesNum = numel(theSourceRGCinputConeIndices);
            if (allRGCsConeInputsNum(theSourceRGCindex) ~= theSourceRGCinputConesNum)
                error('mismatch here')
            end
    
            theSourceRGCinputConeWeights = full(obj.coneConnectivityMatrix(theSourceRGCinputConeIndices,theSourceRGCindex));
            
            % Must be at least 2 cone inputs
            if (theSourceRGCinputConesNum<2)
                continue;
            end
        

            % Source RGC has at least 2 input cones. We will transfer half to
            % the destinationRGC
            destinationRGCindex = zeroInputRGCindices(1);

            % Drop used zero input RGC index
            zeroInputRGCindices = zeroInputRGCindices(2:end);
    
            obj.optimizeTransferOfConeInputsToZeroInputRGC(...
                 theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, destinationRGCindex);
                        
            if (~isempty(videoOBJ))
                % Visualize current connectivity
                hFig = obj.visualizeCurrentConnectivityState(1008);
                videoOBJ.writeVideo(getframe(hFig));
            end
        end % (~isempty(zeroInputRGCindices))

    end % iRGC

    % Update the local RGCRFspacings
    obj.updateLocalRGCRFspacingsBasedOnCurrentCentroids();
end