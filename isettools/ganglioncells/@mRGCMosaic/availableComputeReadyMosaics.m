function availableComputeReadyMosaics(rgcMosaicType)
    validMRGCMosaicTypes = {'ONcenter'};
    assert(ismember(rgcMosaicType, validMRGCMosaicTypes), sprintf('Unknown mRGCmosaic type: ''%s''.', rgcMosaicType));

    p = getpref('isetbio');
    switch (p.rgcResources.method)
       case 'localFile'
           allRGCmosaicsRootDir = p.rgcResources.URLpath;
       otherwise
          error('Unknown rgcResourcesMethod: ''%''.', p.rgcResources.method)
    end % switch
    
    midgetRGCmosaicsRootDir = fullfile(allRGCmosaicsRootDir, sprintf('%sMidgetRGCmosaics', rgcMosaicType),'computeReadyMosaics');
    filesFound = dir(midgetRGCmosaicsRootDir);

    mosaicsAvailableNum = 0;
    for i = 1:numel(filesFound)
        theFileName = filesFound(i).name;
        theFileDir = filesFound(i).folder;
        theFileDate = filesFound(i).date;
        theFileSize = filesFound(i).bytes;
        if (contains(theFileName, 'ComputeReadyMosaic'))
            mosaicsAvailableNum  = mosaicsAvailableNum  + 1;
            theFileDescriptors{mosaicsAvailableNum } = ...
                sprintf('[%02d] %150s. This mosaic was generated on %s and it is %2.1f MBytes.', mosaicsAvailableNum, theFileName, theFileDate, theFileSize/(1024*1024));
        end
    end

    rgcResources = p.rgcResources;
    fprintf('\nThe ''rgcResources'' ISETBio preference is:');
    rgcResources
    
    
    if (mosaicsAvailableNum == 0)
        fprintf(2,'\nDid not find any compute-ready mRGCmosaics of type ''%s'' at %s.\n', rgcMosaicType, midgetRGCmosaicsRootDir);
    else
        fprintf('\nFound %d compute-ready mRGCmosaics of type ''%s'' at %s.\n', mosaicsAvailableNum, rgcMosaicType, midgetRGCmosaicsRootDir);
        fprintf('These are:\n')
        for i = 1:mosaicsAvailableNum
            fprintf('%s\n', theFileDescriptors{i});
        end

end