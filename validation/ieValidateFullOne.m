function ieValidateFullOne(varargin)
% Full data check (no figures, no publish) of pne validation functions
%
%    validateFullOne(param,val, ...)
%
% Possible parameters are:
%    'verbosity' -    high, med, low ...
%    'numeric tolerance' - val
%    'graph mismatched data' - true/false
%    'generate plots' - true, false
%    'doFullAndFastValidation'  - true, false
%
% Examples:
%   validateFullOne('verbosity','high');
%   validateFullOne('Numeric Tolerance',1000*eps);
%   validateFullOne('generate plots',true);
%
% NC, ISETBIO Team, Copyright 2015

%% Close all figures so that we start with a clean slate
close all;

%% We will use preferences for the 'isetbio' project
thisProject = 'isetbio';
UnitTest.usePreferencesForProject(thisProject, 'reset');

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
if (~ispref(thisProject, 'numericTolerance'))
    UnitTest.setPref('numericTolerance', 500*eps);
end

%% Whether to plot data that do not agree with the ground truth
UnitTest.setPref('graphMismatchedData', true);

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
                case 'graphmismatcheddata'
                    UnitTest.setPref('graphMismatchedData', val);
                case 'generateplots'
                    UnitTest.setPref('generatePlots',  val);
                case 'dofullandfastvalidation'
                    if val
                        fullValidationMode = 'FULL';
                    end
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

%% Print all existing validation scripts and ask the user to select one for validation
singleScriptToValidate = UnitTest.selectScriptFromExistingOnes();

%% Validate
UnitTest.runValidationSession({{singleScriptToValidate, []}}, fullValidationMode);

end