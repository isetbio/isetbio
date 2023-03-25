function visualizeVisualRFmapsForMultipleTargetRGCs(...
    theComputeReadyMRGCmosaic, ...
    optimallyMappedSubspaceRFmapsFileName, ...
    mRGCMosaicSubspaceResponsesFileName, ...
    pdfFileName)


    visualizeAnotherSingleRGC = true;
    while (visualizeAnotherSingleRGC)
        
        targetRGCposition = input('Enter (xy) position of target RGC (e.g., [5.6 -1.3]): ');
        targetCenterConesNum = input('Enter # of center cones num (e.g, 3): ');
        targetCenterConeMajorityType = input('Enter type of majority center cone num (either cMosaic.LCONE_ID or cMosaic.MCONE_ID): ');
    
        MosaicPoolingOptimizer.visualizeVisualRFmapForTargetRGC(...
                theComputeReadyMRGCmosaic, ...
                optimallyMappedSubspaceRFmapsFileName, ...
                targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
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
    