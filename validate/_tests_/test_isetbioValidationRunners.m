function tests = test_isetbioValidationRunners()
tests = functiontests(localfunctions);
end

function testPublicRunnerLocations(testCase)
%% Public validation entry points live together in validate.

validateDir = fileparts(fileparts(mfilename('fullpath')));
runnerNames = {'isetbioUnitTest','isetbioTutorialTest','isetbioExampleTest'};
for runnerIndex = 1:numel(runnerNames)
    runnerPath = which(runnerNames{runnerIndex});
    verifyEqual(testCase,fileparts(runnerPath),validateDir);
end

end

function testTutorialOptionNames(testCase)
%% Removed option aliases fail before a tutorial run begins.

verifyError(testCase,@() isetbioTutorialTest('select','t_missing'), ...
    'isetbioTutorialTest:InvalidOption');
verifyError(testCase,@() isetbioTutorialTest('t_missing'), ...
    'isetbioTutorialTest:InvalidInput');

end

function testExampleOptionNames(testCase)
%% Removed option aliases fail before an example run begins.

verifyError(testCase,@() isetbioExampleTest('select','s_missing'), ...
    'isetbioExampleTest:InvalidOption');
verifyError(testCase,@() isetbioExampleTest('s_missing'), ...
    'isetbioExampleTest:InvalidInput');

end
