function generateLocalPrefs()

    resourcesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', fullfile(resourcesRootDir, 'IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM'), ...
        'rawFigurePDFsDir',    fullfile(resourcesRootDir, 'ManuscriptSupportMaterials/PLOS2024/figures/SLIM/raw') ...
        );
     
     % Set the rgcResources pref
     setpref('isetbio', 'rgcResources', rgcResources);

end
