function validateDemo
%
% Validation demo illustrating how to 
% - set various validation preferences
% - validate a list of scripts or a list of script directories. 
% - conduct a validationSession with different modes
    
    %% We will use preferences for the 'isetbioValidation' project - this is project specific
    UnitTest.usePreferencesForProject('isetbioValidation', 'reset');

    %% Change some preferences:
    %% Run time error behavior
    % valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
    UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');
    
    %% Plot generation
    UnitTest.setPref('generatePlots',  false); 
    UnitTest.setPref('closeFigsOnInit', true);
    
    %% Verbosity Level
    % valid options are: 'none', min', 'low', 'med', 'high', 'max'
    UnitTest.setPref('verbosity', 'med');
    
    %% Numeric tolerance for comparison to ground truth data
    UnitTest.setPref('numericTolerance', 500*eps);
    
    %% Whether to plot data that do not agree with the ground truth
    UnitTest.setPref('graphMismatchedData', true);
    
    %% Print current values of isetbioValidation prefs
    UnitTest.listPrefs();
    
    %% What to validate
    validateAllDirs = false;
    if (validateAllDirs)
        % List of script directories to validate. Each entry contains a cell array with 
        % with a validation script directory and an optional struct with
        % prefs that override the corresponding isetbioValidation prefs.
        % At the moment only the 'generatePlots' pref can be overriden.
        %
        % We have functions that generate various stock vScriptsList, but
        % here we do it out explicitly.
        % 
        % Note that the exampleScripts directory contains some scripts that
        % intentionally fail in various ways.
        
        % Get rootDir
        rootDir = UnitTest.getPref('validationRootDir');

        vScriptsList = {...
                {fullfile(rootDir, 'scripts', 'color')} ... 
                {fullfile(rootDir, 'scripts', 'cones')} ... 
                {fullfile(rootDir, 'scripts', 'human')} ... 
                {fullfile(rootDir, 'scripts', 'optics')} ... 
                {fullfile(rootDir, 'scripts', 'radiometry')} ... 
                {fullfile(rootDir, 'scripts', 'scene')} ... 
                {fullfile(rootDir, 'scripts', 'codedevscripts'), struct('generatePlots', true) } ...
            };
    else
        % Alternatively, you can provide a list of scripts to validate. 
        % In this case each entry contains a cell array with 
        % with a script name and an optional struct with
        % prefs that override the corresponding isetbioValidation prefs.
        % At the moment only the generatePlots pref can be overriden.
        
        % Get rootDir
        rootDir = UnitTest.getPref('validationRootDir');
        
        vScriptsList = {...
           {fullfile(rootDir, 'scripts', 'color', 'v_stockman2xyz.m')}
           {fullfile(rootDir, 'scripts', 'codedevscripts', 'v_skeleton.m'), struct('generatePlots', true) }
        };
    end
    
    %% How to validate
    % Run a RUN_TIME_ERRORS_ONLY validation session
    % UnitTest.runValidationSession(vScriptsList, 'RUN_TIME_ERRORS_ONLY')
    
    % Run a FAST validation session (comparing SHA-256 hash keys of the data)
    %UnitTest.runValidationSession(vScriptsList, 'FAST');
    
    % Run a FULL validation session (comparing actual data)
    % UnitTest.runValidationSession(vScriptsList, 'FULL');
    
    % Run a PUBLISH validation session (comparing actual data and update github wiki)
    % UnitTest.runValidationSession(vScriptsList, 'PUBLISH);
    
    % Run a validation session without a specified mode. You will be
    % promped to select one of the available modes.
    UnitTest.runValidationSession(vScriptsList);
end