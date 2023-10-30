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
%    You should be able just to use this template version without editing
%    for both full and fast validations -- it pulls its data off a server
%    we have set up via the RemoteDataToolbox facility.
%
%    Note that you still want to have the remoteDataToolboxConfig set up
%    even if you aren't using the RemoteDataToolbox for validations.
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
% 'useRemoteDataToolbox'    - Boolean. If true use the Remote Data Toolbox
%                             to fetch full validation data on demand.
%                             Default false.
% 'remoteDataToolboxConfig' - String. The struct, file path, or project
%                             name with Remote Data Toolbox configuration.
%                             Default 'isetbio'.
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
    otherwise
        % Some unspecified machine, try user specific customization
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

        % switch(computerInfo.userShortName)
        %     % Could put user specific things in, but at the moment generic
        %     % is good enough.
        % 
        %     case 'colorlab'
        %         % SACCSFA desktop (Linux)
        %         userNameDropbox = 'Mela Nopsin';
        %         dropboxValidationRootDirPath = fullfile('/home/',computerInfo.userShortName,'Aguirre-Brainard Lab Dropbox',userNameDropbox);
        % 
        %     otherwise
        %         dropboxValidationRootDirPath = fullfile('/Users/',computerInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
        % end
end

% RGC mosaic resources Dropbox URLpath
if (~isempty(dropboxValidationRootDirPath))
    rgcDropboxURLpath = fullfile(dropboxValidationRootDirPath, 'IBIO_rgcMosaicResources');
else
    rgcDropboxURLpath = '';
end

if (exist('isetvalidateRootPath','file'))
    validationRootDir = fullfile(isetvalidateRootPath, 'isetbio');
else
    validationRootDir = fullfile(isetbioRootPath, 'validation');
end

p = struct(...
    'projectName', 'isetbio', ...
    'validationRootDir', validationRootDir, ...
    'alternateFastDataDir', '', ...
    'alternateFullDataDir', fullfile(dropboxValidationRootDirPath, 'ISETBioValidationFiles/gradleFiles/validationFull'), ...
    'useRemoteDataToolbox', ~true, ...
    'remoteDataToolboxConfig', 'isetbio', ...
    'githubRepoURL', 'http://isetbio.github.io/isetbio', ...
    'rgcResources', struct('method', 'localFile', 'URLpath', rgcDropboxURLpath), ...
    'generateGroundTruthDataIfNotFound', ~true, ...
    'listingScript', 'ieValidateListAllValidationDirs', ...
    'coreListingScript', 'ieValidateListCoreValidationFiles', ...
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