function isetbioLocalHook
% isetbioLocalHook
%
% Method to set ISETBIO-specific preferences. 
%
% Generally, this function should be edited for your site and then run once.
%
% You should be able just to use this template version without editing
% for both full and fast validations -- it pulls its data off a server we
% have set up via the RemoteDataToolbox facility.
%
% Note that you still want to have the remoteDataToolboxConfig set up even
% if you aren't using the RemoteDataToolbox for validations.
%
% If you are using the ToolboxToolbox, the copy this file to your local
% hooks folder and delete "Template" from the filename.  The TbTb will then
% run it when you deploy isetbio.

    % Specify project-specific preferences for unit test toolbox
    p = struct(...
            'projectName',           'isetbio', ...                                                                                   % The project's name (also the preferences group name)
            'validationRootDir',     fullfile(isetbioRootPath, 'validation'), ...                                                     % Directory location where the 'scripts' subdirectory resides.
            'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
            'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
            'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
            'remoteDataToolboxConfig', 'isetbio', ...                                                                                 % Struct, file path, or project name with Remote Data Toolbox configuration.
            'githubRepoURL',         'http://isetbio.github.io/isetbio', ...                                                          % Github URL for the project. This is only used for publishing tutorials.
            'generateGroundTruthDataIfNotFound',   false, ...                                                                         % Flag indicating whether to generate ground truth if one is not found
            'listingScript',         'ieValidateListAllValidationDirs', ...    
            'coreListingScript',     'ieValidateListCoreValidationFiles', ...
            'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
        );

    generatePreferenceGroup(p);
    UnitTest.usePreferencesForProject(p.projectName);
end

function generatePreferenceGroup(p)
    % remove any existing preferences for this project
    if ispref(p.projectName)
        rmpref(p.projectName);
    end
    
    % generate and save the project-specific preferences
    setpref(p.projectName, 'projectSpecificPreferences', p);
    fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end