function tests = test_outerSegmentUtilities()
tests = functiontests(localfunctions);
end

function testResampleTwoAndThreeDimensionalSignals(testCase)
originalTime = [0 1 2];
newTime = [0 0.5 1 1.5 2];
signal2D = [0 1 2; 2 1 0];
signal3D = reshape(1:12, [2 2 3]);

resampled2D = outerSegment.resample(signal2D, originalTime, newTime);
resampled3D = outerSegment.resample(signal3D, originalTime, newTime);

testCase.verifySize(resampled2D, [2 5]);
testCase.verifyEqual(resampled2D(:, [1 3 5]), signal2D, 'AbsTol', 1e-12);
testCase.verifySize(resampled3D, [2 2 5]);
testCase.verifyEqual(resampled3D(:, :, [1 3 5]), signal3D, 'AbsTol', 1e-12);
end

function testFrozenNoiseIsRepeatable(testCase)
absorptions = reshape([1 5 25 100], [2 2]);
[first, firstNoise] = coneMosaicRect.photonNoise(absorptions, 'noiseFlag', 'frozen', 'seed', 17);
[second, secondNoise] = coneMosaicRect.photonNoise(absorptions, 'noiseFlag', 'frozen', 'seed', 17);

testCase.verifyEqual(first, second);
testCase.verifyEqual(firstNoise, secondNoise);
testCase.verifyGreaterThanOrEqual(first, zeros(size(first)));
testCase.verifyEqual(first, [0 17; 3 81]);
end
