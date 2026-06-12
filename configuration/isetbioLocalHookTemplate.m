function isetbioLocalHookTemplate 
% Method to set ISETBIO-specific preferences. 
%
% Syntax:
%   isetbioLocalHookTemplate
%
% Description:
%    Generally, this function should be edited for your website and then
%    run once.
%
%    You should be able to use this template version without editing for
%    both full and fast validations when the validation data are available
%    in the configured local directory.
%
%    If you are using the ToolboxToolbox, the copy this file to your local
%    hooks folder and delete "Template" from the filename.  The TbTb will
%    then run it when you deploy isetbio.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Specify project-specific preferences for unit test toolbox
% 'projectName'             - String. The project's name (also the
%                             preferences group name). Default 'isetbio'
% 'validationRootDir'       - String. The directory location where the
%                             'scripts' subdirectory resides. Default
%                             fullfile(isetbioRootPath, 'validation')
% 'alternateFastDataDir'    - String. An alternate FAST (hash) data
%                             directory location. Specify '' to use the
%                             default location, i.e.,
%                             $validationRootDir/data/fast. Default ''.
% 'alternateFullDataDir'    - String. An alternate FULL data directory
%                             location. Specify '' to use the default
%                             location, i.e., $validationRootDir/data/full.
%                             Default ''.
% 'githubRepoURL'           - String. The Github URL for the project. This
%                             is only used for publishing tutorials.
%                             Default 'http://isetbio.github.io/isetbio'.
% 'generateGroundTruthDataIfNotFound'
%                           - Boolean. Flag indicating whether or not to
%                             generate the ground truth if one is not
%                             found. Default false.
% 'listingScript'           - String. The listing script name. Default
%                             'ieValidateListAllValidationDirs'.
% 'coreListingScript'       - String. The core listing script name. Default
%                             'ieValidateListCoreValidationFiles'.
% 'numericTolerance'        - Numeric. Numeric tolerance for comparisons
%                             with validation data. Default 1e-11.


% Get Dropbox Validation RootDir location
computerInfo = GetComputerInfo();

if (strcmp(computerInfo.MatlabPlatform, 'GLNXA64'))
    % In Linux, usr networkName instead of localHostName
    computerInfo.localHostName = computerInfo.networkName;
end

switch (computerInfo.localHostName)
    case 'Leviathan'
        % Leviathan
        dropboxValidationRootDirPath = '/mnt/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
    case 'Ithaka'
        % Nicolas' M1 Macbook Pro
        dropboxValidationRootDirPath = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
    case 'Crete'
        % Nicolas' M1 MacStudio Ultra
        dropboxValidationRootDirPath = '/Volumes/M1ProBackUp/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
    case 'Santorini'
        % Nicolas' M1 MacMini
        dropboxValidationRootDirPath  = '/Users/nicolas/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';    
    otherwise
         if ismac
            dbJsonConfigFile = '~/.dropbox/info.json';
            fid = fopen(dbJsonConfigFile);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            dropboxValidationRootDirPath = val.business.path;
            %dropboxValidationRootDirPath = fullfile('/Users/',sysInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
         else
            error('Dropbox validation root directory location not available for computer named: ''%s''.', computerInfo.localHostName);
         end
end

% RGC mosaic resources Dropbox URLpath
if (~isempty(dropboxValidationRootDirPath))
    rgcDropboxURLpath = fullfile(dropboxValidationRootDirPath, 'IBIO_rgcMosaicResources');
else
    rgcDropboxURLpath = '';
end

if (exist('isetvalidateRootPath','file'))
    % Preserve the external isetvalidate repository's legacy directory and
    % listing-script names until that repository migrates them.
    validationRootDir = fullfile(isetvalidateRootPath, 'isetbioRDT');
    listingScript = 'ieValidateRDTListAllValidationDirs';
else
    validationRootDir = fullfile(isetbioRootPath, 'validation');
    listingScript = 'ieValidateListAllValidationDirs';
end

% UnitTestToolbox currently reads the two disabled RemoteDataToolbox fields
% unconditionally, so retain them as inert compatibility preferences.
p = struct(...
    'projectName', 'isetbio', ...
    'validationRootDir', validationRootDir, ...
    'alternateFastDataDir', '', ...
    'alternateFullDataDir', fullfile(dropboxValidationRootDirPath, 'ISETBioValidationFiles/gradleFiles/validationFull'), ...
    'useRemoteDataToolbox', false, ...
    'remoteDataToolboxConfig', '', ...
    'githubRepoURL', 'http://isetbio.github.io/isetbio', ...
    'rgcResources', struct('method', 'localFile', 'URLpath', rgcDropboxURLpath), ...
    'generateGroundTruthDataIfNotFound', true, ...
    'listingScript', listingScript, ...
    'coreListingScript', '', ...
    'numericTolerance', 1e-11);

% Add to the path the Dropbox Validation RootDir location
addpath(genpath(p.alternateFullDataDir));

generatePreferenceGroup(p);
UnitTest.usePreferencesForProject(p.projectName);

% Some behaviors are controlled through the ISET preferences.
%
% Turn off the waitbar!
setpref('ISET','waitbar',false);
ieSessionSet('waitbar',false);


end

function generatePreferenceGroup(p)
% Remove any existing preferences for this project
%
% Syntax:
%   generatePreferenceGroup(p)
%
% Description:
%    Remove existing preferences for the provided project p.
%
% Inputs:
%    p - Object. A project object.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

if ispref(p.projectName), rmpref(p.projectName); end

% generate and save the project-specific preferences
setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf('Generated and saved preferences specific to the %s project.\n', p.projectName);
end
