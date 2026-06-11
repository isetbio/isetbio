function tests = test_cMosaic()
tests = functiontests(localfunctions);
end

function testOuterSegmentLengthModel(testCase)
[lengths, fovealLength] = cMosaic.outerSegmentLengthFromEccentricity([0 4.86 39.90]);

testCase.verifyEqual(fovealLength, 47.81, 'AbsTol', 1e-12);
testCase.verifyEqual(lengths, [47.81 21.2 13.22], 'AbsTol', 1e-12);
testCase.verifyGreaterThanOrEqual(lengths(1:end-1), lengths(2:end));
end

function testNoisyInstancesFrozenSeed(testCase)
meanAbsorptions = reshape([1 10 30 100], [2 2]);
first = cMosaic.noisyInstances(meanAbsorptions, 'noiseFlag', 'frozen', 'seed', 9);
second = cMosaic.noisyInstances(meanAbsorptions, 'noiseFlag', 'frozen', 'seed', 9);

testCase.verifyEqual(first, second);
testCase.verifyGreaterThanOrEqual(first, zeros(size(first)));
end
