function tests = test_fixationalEM()
tests = functiontests(localfunctions);
end

function testDefaultParameters(testCase)
em = fixationalEM;

testCase.verifyEqual(em.timeStepDurationSeconds, 0.001, 'AbsTol', 1e-12);
testCase.verifyEqual(em.microSaccadeType, 'heatmap/fixation based');
testCase.verifyEqual(em.microSaccadeMeanIntervalSeconds, 0.45, ...
    'AbsTol', 1e-12);
testCase.verifyEqual(em.feedbackXposDelaySeconds, 0.07, 'AbsTol', 1e-12);
testCase.verifyEqual(em.feedbackYposDelaySeconds, 0.04, 'AbsTol', 1e-12);
testCase.verifyFalse(em.displayComputeProgress);
testCase.verifyFalse(em.beVerbose);
end

function testSeededComputeIsRepeatable(testCase)
first = configuredFixationalEM(7);
second = configuredFixationalEM(7);

first.compute(0.02, 0.005, 2, true);
second.compute(0.02, 0.005, 2, true);

testCase.verifyEqual(first.timeAxis, 0:0.005:0.015, 'AbsTol', 1e-12);
testCase.verifySize(first.emPosArcMin, [2 4 2]);
testCase.verifySize(first.velocityArcMin, [2 4]);
testCase.verifyEqual(first.emPosArcMin, second.emPosArcMin, 'AbsTol', 1e-12);
testCase.verifyEqual(first.velocityArcMin, second.velocityArcMin, ...
    'AbsTol', 1e-12);
testCase.verifyGreaterThanOrEqual(first.velocityArcMin, ...
    zeros(size(first.velocityArcMin)));
end

function testPathCentering(testCase)
em = configuredFixationalEM(11);
em.compute(0.02, 0.005, 2, false, 'centerPaths', true);

pathMeans = mean(em.emPosArcMin, 2);
testCase.verifyEqual(pathMeans, zeros(size(pathMeans)), 'AbsTol', 1e-12);
testCase.verifyEmpty(em.velocityArcMin);
end

function testCenterAtSpecificTime(testCase)
em = configuredFixationalEM(13);
em.compute(0.02, 0.005, 2, false, ...
    'centerPathsAtSpecificTimeMsec', 10);

[~, index] = min(abs(1000 * em.timeAxis - 10));
positionsAtCenterTime = em.emPosArcMin(:, index, :);
testCase.verifyEqual(positionsAtCenterTime, ...
    zeros(size(positionsAtCenterTime)), 'AbsTol', 1e-12);
end

function testFixationMap(testCase)
timeAxis = 0:0.1:0.3;
paths = zeros(1, 4, 2);
paths(1, :, 1) = [-0.75 -0.25 0.25 0.75];
paths(1, :, 2) = [-0.75 -0.25 0.25 0.75];

[map, supportX, supportY, xSlice, ySlice] = ...
    fixationalEM.computeFixationMap(timeAxis, paths, [-1 1], 1);

testCase.verifySize(map, [3 3]);
testCase.verifyEqual(max(map(:)), 1, 'AbsTol', 1e-12);
testCase.verifyEqual(supportX, [-0.5 0.5 1.5], 'AbsTol', 1e-12);
testCase.verifyEqual(supportY, supportX, 'AbsTol', 1e-12);
testCase.verifyEqual(xSlice, [0 1 0], 'AbsTol', 1e-12);
testCase.verifyEqual(ySlice, [0; 1; 0], 'AbsTol', 1e-12);
end

function testRectangularMosaicGenerationUsesFixationalEM(testCase)
mosaic = coneMosaicRect('pattern', [2 3; 4 2], ...
    'noiseFlag', 'none', 'useParfor', false);

[path, em] = mosaic.emGenSequence(4, 'nTrials', 2, 'rSeed', 7);

testCase.verifyClass(em, 'fixationalEM');
testCase.verifySize(path, [2 4 2]);
testCase.verifyEqual(mosaic.emPositions, path);
end

function em = configuredFixationalEM(seed)
em = fixationalEM;
em.stabilizationSeconds = 0.1;
em.microSaccadeType = 'none';
em.randomSeed = seed;
end
