function visualizeVisualRFmapsForMultipleTargetRGCs(...
    theComputeReadyMRGCmosaic, ...
    optimallyMappedSubspaceRFmapsFileName, ...
    mRGCMosaicSubspaceResponsesFileName, ...
    tickSeparationArcMin, spatialSupportRangeArcMin, ...
    pdfFileName)


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


        MosaicPoolingOptimizer.visualizeVisualRFmapForTargetRGC(...
                theComputeReadyMRGCmosaic, ...
                optimallyMappedSubspaceRFmapsFileName, ...
                targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
                tickSeparationArcMin, spatialSupportRangeArcMin, ...
                pdfFileName);

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
            mRGCMosaicSubspaceResponsesFileName);
    end
    