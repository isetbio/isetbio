function performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsToInclude)

   employTemporalEquivalentEccentricity = true;

   theMosaicFileNames = cell(1, numel(mosaicEccsToInclude));
   theMRGCSTFResponsesFileNames = cell(1, numel(mosaicEccsToInclude));

   % Generate all filenames for all included mosaics
   for iMosaic = 1:numel(mosaicEccsToInclude)
       mosaicParams = MosaicPoolingOptimizer.getMosaicParams(mosaicEccsToInclude(iMosaic));

       % Ask the user which optics were used for computing the input cone
       % mosaic STF responses, so we can obtain the corresponding coneMosaicSTFresponsesFileName
       opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

       % Ask the user which H1 cell index to use for optimizing the RF
       % surround pooling model
       retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

       % Generate the filename of the compute-ready mRGCMosaic
       [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

      % Generate filename for the computed mRGCMosaicSTF responses
      [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

      theMosaicFileNames{iMosaic} = fullfile(resourcesDirectory, computeReadyMosaicFileName);
      theMRGCSTFResponsesFileNames{iMosaic} = fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName);
   end % iEcc

   MosaicPoolingOptimizer.visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
            theMosaicFileNames, ...
            theMRGCSTFResponsesFileNames, ...
            employTemporalEquivalentEccentricity);

end

