function validateFastAll
%
% Fast run (no figures, hash data check, no publish) of all of our
% validation functions.

    close all
    clc
    
    %% We will use preferences for the 'isetbioValidation' project - this is project specific
    UnitTest.usePreferencesForProject('isetbioValidation', 'reset');                
                                       
    %% Set default preferences for this function
    
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
    UnitTest.setPref('graphMismatchedData', false);

    %% Print current values of isetbioValidation prefs
    UnitTest.listPrefs();
    
    %% What to validate
    vScriptsList = validateListAllValidationDirs;
    
    %% How to validate
    % Run a FAST validation session (comparing SHA-256 hash keys of the data)
    UnitTest.runValidationSession(vScriptsList, 'FAST');

end