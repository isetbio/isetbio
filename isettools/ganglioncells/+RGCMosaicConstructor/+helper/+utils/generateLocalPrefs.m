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
            intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics');
            figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/denovo/PLOS2024/figures');

        case 'Crete'
            % User on Crete (Nicolas' M2 Mac Studio) 
            rgcResourcesRootDir  = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
            intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics');
            figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/denovo/PLOS2024/figures');

        % Different computer names for different users here
        case 'YOUR_COMPUTER_NAME'
            %rgcResourcesRootDir  = '...';
            %intermediateDataDir = fullfile(rgcResourcesRootDir,'...');
            %figurePDFsDir = fullfile(rgcResourcesRootDir, '...');

        otherwise
            error('No rgcResourcesRootDir specified for computer with name: ''%s''. Edit generateLocalPrefs.m to generate specific ', p.localHostName);
    end


    % User needs to make sure that the following dirs exist under the rgcResourcesRootDir
    if (useSLIMpaths)
        intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM');
        figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/PLOS2024/figures/SLIM/raw');
    else
        
    end

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', intermediateDataDir, ...
        'figurePDFsDir',    figurePDFsDir ...
     );
     
     % Set the rgcResources pref
     setpref('isetbio', 'rgcResources', rgcResources);
end
