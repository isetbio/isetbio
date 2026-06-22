function run = isetbioTutorialsTest(varargin)
% Run ISETBio tutorials through the shared ISETCam test engine.
%
% Syntax:
%   run = isetbioTutorialsTest
%   run = isetbioTutorialsTest('select',scriptName)
%   run = isetbioTutorialsTest('start',scriptName)
%
% With no arguments, all tutorials run.  'select' runs only scriptName;
% 'start' runs scriptName and every tutorial after it.  A single scriptName
% argument remains supported as the legacy form of 'select'.

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
%% Parse the public selection options while retaining legacy calls.

selector = '';
start = '';
if isempty(varargin), return; end
if isscalar(varargin)
    selector = varargin{1};
    return;
end
if numel(varargin) ~= 2
    error('isetbioTutorialsTest:InvalidInput', ...
        'Use no arguments or one name-value pair: select or start.');
end

option = lower(char(varargin{1}));
switch option
    case {'select','selector'}
        selector = varargin{2};
    case 'start'
        start = varargin{2};
    otherwise
        error('isetbioTutorialsTest:InvalidOption', ...
            'Unknown option "%s". Use select or start.',option);
end

end

function localEnsureISETCam(repoRoot)
%% Add the sibling ISETCam dependency when the shared engine is unavailable.

if ~isempty(which('ieRunTutorialExampleTests')), return; end
dependencyRoot = fullfile(fileparts(repoRoot),'isetcam');
if ~isfolder(dependencyRoot)
    error('isetbioTutorialsTest:MissingISETCam', ...
        'ISETCam dependency not found: %s',dependencyRoot);
end
addpath(genpath(dependencyRoot));

end
