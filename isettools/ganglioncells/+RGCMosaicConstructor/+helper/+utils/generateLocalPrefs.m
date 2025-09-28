% RGCMosaicConstructor.helper.utils.generateLocalPrefs(['useSLIMpaths', true])
function generateLocalPrefs(varargin)
    ip = inputParser;
    ip.addParameter('useSLIMpaths', false, @islogical);
    % Execute the parser
    ip.parse(varargin{:});
    useSLIMpaths = ip.Results.useSLIMpaths;

    p = GetComputerInfo;
    
    % Retrieve the location of the rgcResourcesRootDir. 
    % This is computer/user specific
    switch (p.localHostName)
        case 'Lefkada'
            % User on Lefkada (Nicolas' laptop)
            rgcResourcesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Crete'
            % User on Crete (Nicolas' M2 Mac Studio) 
            rgcResourcesRootDir  = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        % add different computer names for individual users here
        otherwise
            error('No rgcResourcesRootDir for computer with name: ''%s''.', p.localHostName);
    end


    % User needs to make sure that the following dirs exist under the rgcResourcesRootDir
    if (useSLIMpaths)
        intermediateDataDir = 'IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM';
        figurePDFsDir = 'ManuscriptSupportMaterials/PLOS2024/figures/SLIM/raw';
    else
        intermediateDataDir = 'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics/';
        figurePDFsDir = 'ManuscriptSupportMaterials/denovo/PLOS2024/figures/';
    end

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', fullfile(rgcResourcesRootDir, intermediateDataDir), ...
        'figurePDFsDir',    fullfile(rgcResourcesRootDir, figurePDFsDir) ...
     );
     
     % Set the rgcResources pref
     setpref('isetbio', 'rgcResources', rgcResources);
end
