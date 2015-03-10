function vScriptsList = validateListAllValidationDirs
%
% This encapsulates a vull list of our validation directories, so we only
% need to update it in one place.
% 
% Doesn't list the example scripts, and doesn't override any default prefs.
%
% ISETBIO Team (c) 2014

% List of script directories to validate. Each entry contains a cell array with 
% with a validation script directory and an optional struct with
% prefs that override the corresponding isetbioValidation prefs.
% At the moment only the 'generatePlots' pref can be overriden.
%        

% Get rootDir
rootDir = UnitTest.getPref('validationRootDir');

vScriptsList = {...
        {fullfile(rootDir, 'scripts', 'color')} ... 
        {fullfile(rootDir, 'scripts', 'cones')} ... 
        {fullfile(rootDir, 'scripts', 'human')} ... 
        {fullfile(rootDir, 'scripts', 'optics')} ... 
        {fullfile(rootDir, 'scripts', 'radiometry')} ... 
        {fullfile(rootDir, 'scripts', 'scene')} ... 
    };
    
end