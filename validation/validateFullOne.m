function validateFullOne(varargin)
    
    %% We will use preferences for the 'isetbioValidation' project - this is project specific
    UnitTest.usePreferencesForProject('isetbioValidation', 'reset');

    % Run time error behavior
    % valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
    UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');

    % Plot generation
    UnitTest.setPref('generatePlots',  false);
    UnitTest.setPref('closeFigsOnInit', true);

    %% Verbosity Level
    % valid options are: 'none', min', 'low', 'med', 'high', 'max'
    UnitTest.setPref('verbosity', 'high');

    %% Numeric tolerance for comparison to ground truth data
    UnitTest.setPref('numericTolerance', 500*eps);

    %% Whether to plot data that do not agree with the ground truth
    UnitTest.setPref('graphMismatchedData', true);

    %% Print all existing validation scripts and ask the user to select one for validation
    singleScriptToValidate = UnitTest.selectScriptFromExistingOnes();
    
    %% Validate
    UnitTest.runValidationSession({{singleScriptToValidate, []}}, 'FULL');
    
end