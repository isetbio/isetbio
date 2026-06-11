function results = conesUnitTest(mode)
% CONESUNITTEST - Run core cone pigment tests.

if ieNotDefined('mode'), mode = 'core'; end
results = localRunTests(mode, 'conesUnitTest');

end

function results = localRunTests(mode, runnerName)
mode = ieParamFormat(mode);
[testDir, ~, ~] = fileparts(mfilename('fullpath'));
import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;

existingFigures = findall(groot, 'Type', 'figure');
cleanupFigures = onCleanup(@() localCloseTestFigures(existingFigures));
suite = TestSuite.fromFolder(testDir);

switch mode
    case {'core', 'fast', 'quantitative'}
        names = {suite.Name};
        suite = suite(~contains(names, 'FullOnly'));
    case {'full', 'all'}
    otherwise
        error('Unknown %s mode %s. Use ''core'' or ''full''.', runnerName, mode);
end

runner = TestRunner.withTextOutput;
results = runner.run(suite);
ieTestReport(results, runnerName);
end

function localCloseTestFigures(existingFigures)
allFigures = findall(groot, 'Type', 'figure');
testFigures = setdiff(allFigures, existingFigures);
testFigures = testFigures(ishghandle(testFigures));
if ~isempty(testFigures), close(testFigures); end
end
