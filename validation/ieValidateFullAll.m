function ieValidateFullAll(varargin)
% Full data check (no figures, no publish) of all validation functions
%
%    ieValidateFullAll(param,val, ...)
%
% Possible parameters are:
%    'verbosity' -    high, med, low ...
%    'numeric tolerance' - val
%    'graph mismatched data' - true/false
%    'generate plots' - true, false
%    'doFullAndFastValidation'  - true, false
%    'asAssertion' - true/false
%
% Examples:
%   ieValidateFullAll('verbosity','high');
%   ieValidateFullAll('Numeric Tolerance',1000*eps);
%   ieValidateFullAll('generate plots',true);
%
% NC, ISETBIO Team, Copyright 2015

%% Close all figures so that we start with a clean slate
close all;

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'isetbio';
UnitTest.usePreferencesForProject(thisProject, 'reset');

%% Set preferences for this function

% Run time error behavior
% valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');

% Plot generation
UnitTest.setPref('generatePlots',  false);
UnitTest.setPref('closeFigsOnInit', true);

%% Verbosity Level
% Valid options are: 'none', min', 'low', 'med', 'high', 'max'
UnitTest.setPref('verbosity', 'low');

%% Numeric tolerance for comparison to ground truth data
if (~ispref(thisProject, 'numericTolerance'))
    UnitTest.setPref('numericTolerance', 500*eps);
end

%% Whether to plot data that do not agree with the ground truth
UnitTest.setPref('graphMismatchedData', true);

%% Whether to throw an error when a validation fails.
asAssertion = false;

%% Adjust parameters based on input arguments
fullValidationMode = 'FULLONLY';
if ~isempty(varargin)
    if ~isodd(length(varargin))
        for ii=1:2:length(varargin)
            param = ieParamFormat(varargin{ii});
            val   = varargin{ii+1};
            switch(param)
                case 'verbosity'
                    UnitTest.setPref('verbosity',val);
                case 'numerictolerance'
                    UnitTest.setPref('numericTolerance', val);
                case 'graphMismatchedData'
                    UnitTest.setPref('graphMismatchedData', val);
                case 'generatePlots'
                    UnitTest.setPref('generatePlots',  val);
                case 'dofullandfastvalidation'
                    if val
                        fullValidationMode = 'FULL';
                    end
                case 'asassertion'
                    asAssertion = val;
                otherwise
                    error('Unknown validation string %s\n',varargin{ii+1});
            end
        end
    else
        error('Odd number of arguments, must be param/val pairs');
    end
end

%% Print current values of isetbioValidation prefs
UnitTest.listPrefs();

%% What to validate
listingScript = UnitTest.getPref('listingScript');
vScriptsList = eval(listingScript);

%% How to validate
% Run a FULL validation session (comparing actual data)
obj = UnitTest.runValidationSession(vScriptsList, fullValidationMode);

if asAssertion
    % Assert no failed validations
    summary = [obj.summaryReport{:}];
    success = ~any([summary.fullFailed]);
    assert(success, 'One or more validations failed.');
end

end