% RGCMosaicConstructor.helper.util.generateLocalPrefs()
function generateLocalPrefs()

    % User needs to specify the rgc resourcesRootDir
    % User on Ithaka
    rgcResourcesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

    % User on Crete
    gcResourcesRootDir  = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';


    % User needs to make sure that the following dirs exist under resourcesRootDir
    intermediateDataDir = 'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics/';
    figurePDFsDir = 'ManuscriptSupportMaterials/denovo/PLOS2024/figures/';

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', fullfile(gcResourcesRootDir, intermediateDataDir), ...
        'figurePDFsDir',    fullfile(gcResourcesRootDir, figurePDFsDir) ...
     );
     
     % Set the rgcResources pref
     setpref('isetbio', 'rgcResources', rgcResources);
end
