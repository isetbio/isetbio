function run = isetbioExampleTest(varargin)
% Run ISETBio examples through the shared ISETCam test engine.
%
% Syntax:
%   run = isetbioExampleTest
%   run = isetbioExampleTest('selection',scriptName)
%   run = isetbioExampleTest('start',scriptName)
%
% With no arguments, all examples run. 'selection' runs only scriptName;
% 'start' runs scriptName and every example after it.

[selector,start] = localParseSelection(varargin{:});

repoRoot = isetbioRootPath;
localEnsureISETCam(repoRoot);

config = struct();
config.repositoryName = 'ISETBio';
config.repositoryRoot = repoRoot;
config.suiteKind = 'examples';
config.runnerName = mfilename;
config.selector = selector;
config.start = start;
config.skipPathPatterns = { ...
    [filesep 'data' filesep], ...
    ['scripts' filesep 'image' filesep 'jpegFiles'], ...
    ['scripts' filesep 'optics' filesep 'chromAb']};

run = ieRunTutorialExampleTests(config);

end

function [selector,start] = localParseSelection(varargin)
%% Parse the public selection options.

selector = '';
start = '';
if isempty(varargin), return; end
if numel(varargin) ~= 2
    error('isetbioExampleTest:InvalidInput', ...
        'Use no arguments or one name-value pair: selection or start.');
end

option = lower(char(varargin{1}));
switch option
    case 'selection'
        selector = varargin{2};
    case 'start'
        start = varargin{2};
    otherwise
        error('isetbioExampleTest:InvalidOption', ...
            'Unknown option "%s". Use selection or start.',option);
end

end

function localEnsureISETCam(repoRoot)
%% Add the sibling ISETCam dependency when the shared engine is unavailable.

if ~isempty(which('ieRunTutorialExampleTests')), return; end
dependencyRoot = fullfile(fileparts(repoRoot),'isetcam');
if ~isfolder(dependencyRoot)
    error('isetbioExampleTest:MissingISETCam', ...
        'ISETCam dependency not found: %s',dependencyRoot);
end
addpath(genpath(dependencyRoot));

end
