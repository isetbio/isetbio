function tests = test_oiSequence()
tests = functiontests(localfunctions);
end

function testConstructorStoresProperties(testCase)
oi1 = oiCreate;
oi2 = oiCreate;
modFn    = [0.0 0.5 1.0 0.5 0.0];
timeAxis = (0:4) * 0.002;

seq = oiSequence(oi1, oi2, timeAxis, modFn, 'composition', 'blend');

testCase.verifyEqual(seq.length, numel(modFn));
testCase.verifyEqual(seq.timeAxis, timeAxis, 'AbsTol', 1e-15);
testCase.verifyEqual(seq.modulationFunction, modFn, 'AbsTol', 1e-15);
testCase.verifyEqual(seq.composition, 'blend');
end

function testScalarTimeAxisExpansion(testCase)
% A scalar time axis is treated as delta-t and expanded to a full vector.
oi1 = oiCreate;
oi2 = oiCreate;
modFn = ones(1, 5);
dt    = 0.001;

seq = oiSequence(oi1, oi2, dt, modFn);

testCase.verifyEqual(seq.timeAxis, dt * (0:4), 'AbsTol', 1e-15);
testCase.verifyEqual(seq.length, 5);
end

function testTimeAxisLengthMismatch(testCase)
oi1 = oiCreate;
oi2 = oiCreate;

localVerifyErrorContains(testCase, ...
    @() oiSequence(oi1, oi2, (1:5) * 0.001, ones(1, 10)), ...
    'Time axis does not match');
end

function testSpatialMismatch(testCase)
% OIs whose photon arrays differ in rows/cols should fail the spatial check.
nWave = oiGet(oiCreate, 'nwave');
oi1 = oiCreate;
oi1 = oiSet(oi1, 'photons', ones(5, 5, nWave));
oi2 = oiCreate;
oi2 = oiSet(oi2, 'photons', ones(3, 3, nWave));

localVerifyErrorContains(testCase, ...
    @() oiSequence(oi1, oi2, (0:2) * 0.001, ones(1, 3)), ...
    'Mismatch');
end

function testFrameAtIndexAdd(testCase)
[seq, nWave] = localBuildSequence('add', [0.5 1.0]);

frame1 = seq.frameAtIndex(1);  % 1 + 0.5*3 = 2.5
frame2 = seq.frameAtIndex(2);  % 1 + 1.0*3 = 4.0

testCase.verifyEqual(oiGet(frame1, 'photons'), 2.5 * ones(3, 3, nWave), 'AbsTol', 1e-12);
testCase.verifyEqual(oiGet(frame2, 'photons'), 4.0 * ones(3, 3, nWave), 'AbsTol', 1e-12);
end

function testFrameAtIndexBlend(testCase)
[seq, nWave] = localBuildSequence('blend', [0.5 0.0]);

frame1 = seq.frameAtIndex(1);  % (1-0.5)*1 + 0.5*3 = 2.0
frame2 = seq.frameAtIndex(2);  % (1-0.0)*1 + 0.0*3 = 1.0

testCase.verifyEqual(oiGet(frame1, 'photons'), 2.0 * ones(3, 3, nWave), 'AbsTol', 1e-12);
testCase.verifyEqual(oiGet(frame2, 'photons'), 1.0 * ones(3, 3, nWave), 'AbsTol', 1e-12);
end

function testFrameAtIndexXor(testCase)
[seq, nWave] = localBuildSequence('xor', [0.0 0.5]);

frame1 = seq.frameAtIndex(1);  % w=0 -> fixed = 1.0
frame2 = seq.frameAtIndex(2);  % w=0.5 -> 0.5*modulated = 1.5

testCase.verifyEqual(oiGet(frame1, 'photons'), 1.0 * ones(3, 3, nWave), 'AbsTol', 1e-12);
testCase.verifyEqual(oiGet(frame2, 'photons'), 1.5 * ones(3, 3, nWave), 'AbsTol', 1e-12);
end

function testMaxEyeMovementsBasicCount(testCase)
% 10 frames × 5 ms per frame ÷ 10 ms integration = 5 eye movements.
oi1 = oiCreate;
oi2 = oiCreate;
nFrames = 10;
dt      = 0.005;
seq = oiSequence(oi1, oi2, (0:(nFrames-1)) * dt, ones(1, nFrames));

n = seq.maxEyeMovementsNumGivenIntegrationTime(0.010);

testCase.verifyEqual(n, 5);
end

function testTimeStep(testCase)
oi1 = oiCreate;
oi2 = oiCreate;
dt  = 0.005;
seq = oiSequence(oi1, oi2, (0:4) * dt, ones(1, 5));

testCase.verifyEqual(seq.timeStep(), dt, 'AbsTol', 1e-15);
end

% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------

function [seq, nWave] = localBuildSequence(composition, weights)
% Create a 3×3 oiSequence with fixed photons=1 and modulated photons=3.
% Use an explicit count: oiGet(oiCreate,'nwave') returns 0 for a fresh OI
% that has no photons, making all arithmetic produce empty arrays.
nWave = 3;
oi1 = oiCreate;
oi1 = oiSet(oi1, 'photons', ones(3, 3, nWave));
oi2 = oiCreate;
oi2 = oiSet(oi2, 'photons', 3 * ones(3, 3, nWave));
timeAxis = (0:(numel(weights)-1)) * 0.001;
seq = oiSequence(oi1, oi2, timeAxis, weights, 'composition', composition);
end

function localVerifyErrorContains(testCase, functionHandle, expectedText)
didError = false;
try
    functionHandle();
catch err
    didError = true;
    testCase.verifySubstring(err.message, expectedText);
end
testCase.verifyTrue(didError, ...
    sprintf('Expected an error containing: %s', expectedText));
end
