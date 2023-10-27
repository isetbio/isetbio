function availableComputeReadyMosaics(rgcMosaicType)
    validMRGCMosaicTypes = {'ONcenterMidgetRGC'};
    assert(ismember(rgcMosaicType, validMRGCMosaicTypes), sprintf('Unknown mRGCmosaic type: ''%s''.', rgcMosaicType));

    p = getpref('isetbio');
    switch (p.rgcResources.method)
       case 'localFile'
           allRGCmosaicsRootDir = p.rgcResources.URLpath;
       otherwise
          error('Unknown rgcResourcesMethod: ''%''.', p.rgcResources.method)
    end % switch
    
    midgetRGCmosaicsRootDir = fullfile(allRGCmosaicsRootDir, sprintf('%smosaics', rgcMosaicType),'computeReadyMosaics');
    filesFound = dir(midgetRGCmosaicsRootDir);

    mosaicsAvailableNum = 0;
    theFileDescriptors = {};
    theFileInfos = {};
    for i = 1:numel(filesFound)
        theFileName = filesFound(i).name;
        theFileDir = filesFound(i).folder;
        theFileDate = filesFound(i).date;
        theFileSize = filesFound(i).bytes;
        if (contains(theFileName, 'ComputeReadyMosaic'))
            mosaicsAvailableNum  = mosaicsAvailableNum  + 1;
            theFileDescriptors{mosaicsAvailableNum} = ...
                sprintf('%150s.', theFileName);
            theFileInfos{mosaicsAvailableNum} = ...
                sprintf(' Birthday: %s, Size:%2.1f_MBytes', theFileDate, theFileSize/(1024*1024));
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
            fprintf('%s %30s\n', strrep(theFileDescriptors{i}, ' ', ''), theFileInfos{i});
        end

end
