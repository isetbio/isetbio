function tests = test_oiArbitrarySequence()
tests = functiontests(localfunctions);
end

function testConstructorStoresProperties(testCase)
oiList   = {oiCreate, oiCreate, oiCreate};
timeAxis = (0:2) * 0.005;

seq = oiArbitrarySequence(oiList, timeAxis);

testCase.verifyEqual(seq.length, numel(oiList));
testCase.verifyEqual(seq.timeAxis, timeAxis, 'AbsTol', 1e-15);
end

function testConstructorRequiresSameLengths(testCase)
oiList   = {oiCreate, oiCreate, oiCreate};
timeAxis = (0:4) * 0.005;  % 5 elements vs 3 OIs

localVerifyThrows(testCase, @() oiArbitrarySequence(oiList, timeAxis));
end

function testFrameAtIndexReturnsCorrectOI(testCase)
oi1 = oiCreate;
oi2 = oiCreate;
oi3 = oiCreate;
% Make them distinguishable by name.
oi1 = oiSet(oi1, 'name', 'first');
oi2 = oiSet(oi2, 'name', 'second');
oi3 = oiSet(oi3, 'name', 'third');

seq = oiArbitrarySequence({oi1, oi2, oi3}, (0:2) * 0.005);

testCase.verifyEqual(oiGet(seq.frameAtIndex(1), 'name'), 'first');
testCase.verifyEqual(oiGet(seq.frameAtIndex(2), 'name'), 'second');
testCase.verifyEqual(oiGet(seq.frameAtIndex(3), 'name'), 'third');
end

function testLengthMatchesTimeAxis(testCase)
nFrames  = 7;
oiList   = arrayfun(@(~) oiCreate, 1:nFrames, 'UniformOutput', false);
timeAxis = (0:(nFrames-1)) * 0.003;

seq = oiArbitrarySequence(oiList, timeAxis);

testCase.verifyEqual(seq.length, nFrames);
end

function testTimeStep(testCase)
oiList   = {oiCreate, oiCreate, oiCreate};
dt       = 0.004;
timeAxis = (0:2) * dt;

seq = oiArbitrarySequence(oiList, timeAxis);

testCase.verifyEqual(seq.timeStep(), dt, 'AbsTol', 1e-15);
end

function testMaxEyeMovementsBasicCount(testCase)
% 10 frames × 5 ms per frame ÷ 10 ms integration = 5 eye movements.
nFrames  = 10;
dt       = 0.005;
oiList   = arrayfun(@(~) oiCreate, 1:nFrames, 'UniformOutput', false);
timeAxis = (0:(nFrames-1)) * dt;

seq = oiArbitrarySequence(oiList, timeAxis);
n   = seq.maxEyeMovementsNumGivenIntegrationTime(0.010);

testCase.verifyEqual(n, 5);
end

% -------------------------------------------------------------------------
% Helper
% -------------------------------------------------------------------------

function localVerifyThrows(testCase, functionHandle)
didError = false;
try
    functionHandle();
catch
    didError = true;
end
testCase.verifyTrue(didError, 'Expected constructor to throw an error.');
end
