function availableComputeReadyMosaics(rgcMosaicType)
    validMRGCMosaicTypes = {'ONcenterMidgetRGC'};
    assert(ismember(rgcMosaicType, validMRGCMosaicTypes), sprintf('Unknown mRGCmosaic type: ''%s''.', rgcMosaicType));

    p = getpref('isetbio');
    if (isempty(p))
        error('Did not find ''isetbio'' preferences');
    end
    assert(~isempty(p), 'Did not find ''isetbio'' preferences');
    assert(isfield(p, 'rgcResources'), 'Did not find ''rgcResources'' field in isetbio prefs\n');
    assert(isfield(p.rgcResources, 'method'), 'Did not find ''method'' field in isetbio prefs.rgcResources\n');

    switch (p.rgcResources.method)
       case 'localFile'
           allRGCmosaicsRootDir = p.rgcResources.URLpath;
       otherwise
          error('Unknown rgcResourcesMethod: ''%''.', p.rgcResources.method)
    end % switch
    
    midgetRGCmosaicsRootDir = fullfile(allRGCmosaicsRootDir, sprintf('%smosaics', rgcMosaicType),'computeReadyMosaics');
    if (~isfolder(midgetRGCmosaicsRootDir))
        fprintf('\n\nThe %s directory with all compute-ready RGC mosaics was not found.\n', midgetRGCmosaicsRootDir);
        if (strfind(midgetRGCmosaicsRootDir, 'Dropbox'))
            fprintf('This directory exists at https://www.dropbox.com/scl/fo/vam7x2b5uwpavh4qj9dm3/h?rlkey=cm9ivsvrphdu67bn77swlb88v&dl=0\n\n');
        end
        error('Unable to find directory: ''%s''.', midgetRGCmosaicsRootDir);
    end


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
