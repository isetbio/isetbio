function run = isetbioTutorialsTest(selector)
% Run ISETBio tutorials through the shared ISETCam test engine.
%
% Syntax:
%   run = isetbioTutorialsTest
%   run = isetbioTutorialsTest(selector)

if nargin < 1, selector = ''; end

repoRoot = isetbioRootPath;
localEnsureISETCam(repoRoot);

config = struct();
config.repositoryName = 'ISETBio';
config.repositoryRoot = repoRoot;
config.suiteKind = 'tutorials';
config.runnerName = mfilename;
config.selector = selector;
config.skipPathPatterns = { ...
    [filesep 'data' filesep], ...
    ['hyperspectral' filesep 'support']};

run = ieRunTutorialExampleTests(config);

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
