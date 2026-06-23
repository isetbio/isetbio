function run = isetbioTutorialTest(varargin)
% Run ISETBio tutorials through the shared ISETCam test engine.
%
% Syntax:
%   run = isetbioTutorialTest
%   run = isetbioTutorialTest('selection',scriptName)
%   run = isetbioTutorialTest('start',scriptName)
%
% With no arguments, all tutorials run. 'selection' runs only scriptName;
% 'start' runs scriptName and every tutorial after it.

[selector,start] = localParseSelection(varargin{:});

repoRoot = isetbioRootPath;
localEnsureISETCam(repoRoot);

config = struct();
config.repositoryName = 'ISETBio';
config.repositoryRoot = repoRoot;
config.suiteKind = 'tutorials';
config.runnerName = mfilename;
config.selector = selector;
config.start = start;
config.skipPathPatterns = { ...
    [filesep 'data' filesep], ...
    ['hyperspectral' filesep 'support']};

run = ieRunTutorialExampleTests(config);

end

function [selector,start] = localParseSelection(varargin)
%% Parse the public selection options.

selector = '';
start = '';
if isempty(varargin), return; end
if numel(varargin) ~= 2
    error('isetbioTutorialTest:InvalidInput', ...
        'Use no arguments or one name-value pair: selection or start.');
end

option = lower(char(varargin{1}));
switch option
    case 'selection'
        selector = varargin{2};
    case 'start'
        start = varargin{2};
    otherwise
        error('isetbioTutorialTest:InvalidOption', ...
            'Unknown option "%s". Use selection or start.',option);
end

end

function localEnsureISETCam(repoRoot)
%% Add the sibling ISETCam dependency when the shared engine is unavailable.

if ~isempty(which('ieRunTutorialExampleTests')), return; end
dependencyRoot = fullfile(fileparts(repoRoot),'isetcam');
if ~isfolder(dependencyRoot)
    error('isetbioTutorialTest:MissingISETCam', ...
        'ISETCam dependency not found: %s',dependencyRoot);
end
addpath(genpath(dependencyRoot));

end
