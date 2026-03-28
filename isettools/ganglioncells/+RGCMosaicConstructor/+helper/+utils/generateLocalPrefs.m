% RGCMosaicConstructor.helper.utils.generateLocalPrefs()
% 
% If your machine is not listed in the switch case below, call this function
% with key-value pairs that specify the full paths to local
% directories in your machine where you want the interemediate data and
% figure PDFs to be saved, e.g.:
% 
%  RGCMosaicConstructor.helper.utils.generateLocalPrefs(...
%     'intermediateDataDir', '/Users/yourName/Desktop/tmpISETBioData', ...
%     'figurePDFsDir', '/Users/yourName/Desktop/tmpISETBioFigures');
%
% See https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cell-(RGC)-mosaics
% section titled ''Configuring ISETBio to access RGC resources and run RGC
% simulations'' for more info.


function generateLocalPrefs(varargin)
    ip = inputParser;
    ip.addParameter('useSLIMpaths', false, @islogical); % Only for OLD mosaics
    ip.addParameter('intermediateDataDir', 'FullPathToLocalDirectoryWhereGeneratedIntermediateDataWillBeSaved'); 
    ip.addParameter('figurePDFsDir', 'FullPathToLocalDirectoryWhereGeneratedFiguresWillBeSaved'); 
    ip.addParameter('queryUserBeforeGeneratingMissingDir', false);

    % Execute the parser
    ip.parse(varargin{:});
    useSLIMpaths = ip.Results.useSLIMpaths;

    p = GetComputerInfo;
    
    switch (p.MatlabPlatform)
        case 'GLNXA64'
            % Linux
            localHostName = p.networkName;
        otherwise
            localHostName = p.localHostName;
    end % switch

    % Retrieve the location of the rgcResourcesRootDir. 
    % This is computer/user specific.
    % New users need to add a case for their computer
    switch (localHostName)
        case 'Lefkada'
            % User on Lefkada (Nicolas' laptop)
            rgcResourcesRootDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

            % Directory where all the mosaic generation intermediate data will be stored
            intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics');
            
            % Directory where all generated figures will be saved
            figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/denovo/PLOS2024/figures');
            
            % I want the code to stop and ask me if it is OK to generate a missing directory
            queryUserBeforeGeneratingMissingDir = true;

        case 'Crete'
            % User on Crete (Nicolas' M2 Mac Studio) 
            rgcResourcesRootDir  = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
            
            % Directory where all the mosaic generation intermediate data will be stored
            intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics');
            
            % Directory where all generated figures will be saved
            figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/denovo/PLOS2024/figures');

            % I want the code to stop and ask me if it is OK to generate a missing directory
            queryUserBeforeGeneratingMissingDir = false;


        case 'Leviathan'
            % User on Crete (Nicolas' Linux box) 
            rgcResourcesRootDir  = localDropboxDir;
            
            % Directory where all the mosaic generation intermediate data will be stored
            intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/denovo/intermediateFiles/ONcenterMidgetRGCmosaics');
            
            % Directory where all generated figures will be saved
            figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/denovo/PLOS2024/figures');

            % I want the code to stop and ask me if it is OK to generate a missing directory
            queryUserBeforeGeneratingMissingDir = false;

        otherwise
            intermediateDataDir = ip.Results.intermediateDataDir;
            figurePDFsDir = ip.Results.figurePDFsDir;
            queryUserBeforeGeneratingMissingDir = ip.Results.queryUserBeforeGeneratingMissingDir;

    end


    % OLD DATA
    if (useSLIMpaths)
        rgcResourcesRootDir = localDropboxDir;
        intermediateDataDir = fullfile(rgcResourcesRootDir,'IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM');
        figurePDFsDir = fullfile(rgcResourcesRootDir,'ManuscriptSupportMaterials/PLOS2024/figures/SLIM/raw');
        queryUserBeforeGeneratingMissingDir = true;
    end

    % Generate the rgcResources struct
    rgcResources = struct(...
        'method', 'localFile', ...
        'intermediateDataDir', intermediateDataDir, ...
        'figurePDFsDir', figurePDFsDir, ...
        'queryUserBeforeGeneratingMissingDir', queryUserBeforeGeneratingMissingDir ...
     );
     
    
    % Set the rgcResources pref
    setpref('isetbio', 'rgcResources', rgcResources);
end
