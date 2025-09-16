% RGCMosaicConstructor.helper.util.generateLocalPrefs()
function generateLocalPrefs()

    % Ithaka
    resourcesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

    % Crete
    resourcesRootDir = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';


    % User needs to make sure that the following dirs exist under resourcesRootDir
    intermediateDataDir = 'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics/';
    figurePDFsDir = 'ManuscriptSupportMaterials/denovo/PLOS2024/figures/';

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', fullfile(resourcesRootDir, intermediateDataDir), ...
        'figurePDFsDir',    fullfile(resourcesRootDir, figurePDFsDir) ...
     );
     
     % Set the rgcResources pref
     setpref('isetbio', 'rgcResources', rgcResources);
end
