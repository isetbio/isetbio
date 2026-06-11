function visualizeVisualRFmapsForMultipleTargetRGCs(...
    theComputeReadyMRGCmosaic, ...
    optimallyMappedSubspaceRFmapsFileName, ...
    pdfFileName, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);

    p.parse(varargin{:});
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;


    visualizeAnotherSingleRGC = true;
    while (visualizeAnotherSingleRGC)
        
        % Ask the user which RGC to look for:
        % position, # of center cones, majority cone type
        rgcSpecification = input('Plot RGC with specific index (1), or RGC at a target position (2) ? ');
        if (rgcSpecification == 1)
            targetRGCposition = [];
            targetCenterConesNum = [];
            targetCenterConeMajorityType = [];
        else
            targetRGCposition = input('Enter (xy) position of target RGC (e.g., [5.6 -1.3]): ');
            targetCenterConesNum = input('Enter # of center cones num (e.g, 3): ');
            targetCenterConeMajorityType = input('Enter type of majority center cone num (either cMosaic.LCONE_ID or cMosaic.MCONE_ID): ');
        end


        MosaicPoolingOptimizer.retrieveAndVisualizeSubspaceRFmapForTargetRGC(...
                theComputeReadyMRGCmosaic, ...
                optimallyMappedSubspaceRFmapsFileName, ...
                targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
                pdfFileName, ...
                'tickSeparationArcMin', tickSeparationArcMin, ...
                'reverseXDir', reverseXDir, ...
                'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions)

        visualizeSingleRGC = input('Visualize visual RF and profiles for another RGC ? [y=YES] : ', 's');
        if (strcmpi(visualizeSingleRGC, 'y'))
            visualizeAnotherSingleRGC = true;
        else
            visualizeAnotherSingleRGC = false;
        end
    end

    visualizeVisualRFcentersOfAllRGCs = input('Visualize the visual RF centers for all RGCs in the mosaic ? [y=YES] : ', 's');
    if (strcmpi(visualizeVisualRFcentersOfAllRGCs , 'y'))
        MosaicPoolingOptimizer.visualizeVisualRFcentersOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName);
    end
    