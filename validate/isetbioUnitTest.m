function results = isetbioUnitTest(mode)
% ISETBIOUNITTEST - Master runner for all ISETBio unit tests
%
% Usage:
%   results = isetbioUnitTest;
%   results = isetbioUnitTest('full');

if ieNotDefined('mode'), mode = 'core'; end
mode = ieParamFormat(mode);

rootPath = isetbioRootPath;
testDirs = dir(fullfile(rootPath, '**', '_tests_'));

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;

existingFigures = findall(groot, 'Type', 'figure');
cleanupFigures = onCleanup(@() localCloseTestFigures(existingFigures));

masterSuite = [];
for ii = 1:numel(testDirs)
    if testDirs(ii).isdir
        folderPath = fullfile(testDirs(ii).folder, testDirs(ii).name);
        folderSuite = TestSuite.fromFolder(folderPath);
        masterSuite = [masterSuite, folderSuite]; %#ok<AGROW>
    end
end

masterSuite = localSelectMode(masterSuite, mode);

if isempty(masterSuite)
    fprintf('No ISETBio tests found for mode ''%s''.\n', mode);
    results = [];
    ieTestReport(results, 'isetbioUnitTest');
    return;
end

runner = TestRunner.withTextOutput;
results = runner.run(masterSuite);
ieTestReport(results, 'isetbioUnitTest');

end

function suite = localSelectMode(suite, mode)
%% Select core or full tests using the full-only filename convention.

switch mode
    case {'core', 'fast', 'quantitative'}
        names = {suite.Name};
        suite = suite(~contains(names, 'FullOnly'));
    case {'full', 'all'}
        % Keep the complete suite.
    otherwise
        error('Unknown isetbioUnitTest mode %s. Use ''core'' or ''full''.', mode);
end

end

function localCloseTestFigures(existingFigures)
%% Close figures opened by tests while preserving pre-existing figures.

allFigures = findall(groot, 'Type', 'figure');
testFigures = setdiff(allFigures, existingFigures);
testFigures = testFigures(ishghandle(testFigures));
if ~isempty(testFigures), close(testFigures); end

end
