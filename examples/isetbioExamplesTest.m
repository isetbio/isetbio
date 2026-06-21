function run = isetbioExamplesTest(selector)
% Run ISETBio examples through the shared ISETCam test engine.
%
% Syntax:
%   run = isetbioExamplesTest
%   run = isetbioExamplesTest(selector)

if nargin < 1, selector = ''; end

repoRoot = isetbioRootPath;
localEnsureISETCam(repoRoot);

config = struct();
config.repositoryName = 'ISETBio';
config.repositoryRoot = repoRoot;
config.suiteKind = 'examples';
config.runnerName = mfilename;
config.selector = selector;
config.skipPathPatterns = { ...
    [filesep 'data' filesep], ...
    ['scripts' filesep 'image' filesep 'jpegFiles'], ...
    ['scripts' filesep 'optics' filesep 'chromAb']};

run = ieRunTutorialExampleTests(config);

end

function localEnsureISETCam(repoRoot)
%% Add the sibling ISETCam dependency when the shared engine is unavailable.

if ~isempty(which('ieRunTutorialExampleTests')), return; end
dependencyRoot = fullfile(fileparts(repoRoot),'isetcam');
if ~isfolder(dependencyRoot)
    error('isetbioExamplesTest:MissingISETCam', ...
        'ISETCam dependency not found: %s',dependencyRoot);
end
addpath(genpath(dependencyRoot));

end
