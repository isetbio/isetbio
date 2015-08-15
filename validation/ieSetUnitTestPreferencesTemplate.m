% Method to set ISETBIO-specific unit testing preferences. 
%
% Generally, this function should be edited for your site and then run
% once.
%
% To use the unit testing framework in different projects, copy this file
% to the projects' validation directory and adapt the p-struct according
% to that project's specifics.  You will want to change isetbioValidation
% and isetbioRootPath and so forth.
%
% If you just run the distributed (Template) version to accept the
% defaults, the FAST validation in isetbio will run.  But the FULL forms of
% validation need the validation data. These can be found here:
%
%  http://????? - To be determined.
%
%
% ISETBIO Team, 2015

function setIsetbioUnitTestPreferencesTemplate

    % Specify project-specific preferences
    p = struct(...
            'projectName',           'isetbioValidation', ...                                                                         % The project's name (also the preferences group name)
            'validationRootDir',     fullfile(isetbioRootPath, 'validation'), ...                                                     % Directory location where the 'scripts' subdirectory resides.
            'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
            'alternateFullDataDir',  fullfile(filesep,'Users1', 'Shared', 'Dropbox', 'ISETBIOFullValidationData'), ...                % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
            'clonedWikiLocation',    fullfile(filesep,'Users',  'Shared', 'Matlab', 'Toolboxes', 'ISETBIO_Wiki', 'isetbio.wiki'), ... % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
            'clonedGhPagesLocation', fullfile(filesep,'Users',  'Shared', 'Matlab', 'Toolboxes', 'ISETBIO_GhPages', 'isetbio'), ...   % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
            'githubRepoURL',         'http://isetbio.github.io/isetbio', ...                                                          % Github URL for the project. This is only used for publishing tutorials.
            'generateGroundTruthDataIfNotFound',   false, ...  % Flag indicating whether to generate ground truth if one is not found
            'listingScript',         'ieValidateListAllValidationDirs', ...    
            'coreListingScript',     'ieValidateListCoreValidationFiles' ...
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