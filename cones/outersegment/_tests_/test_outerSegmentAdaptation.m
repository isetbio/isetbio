function tests = test_outerSegmentAdaptation()
tests = functiontests(localfunctions);
end

function testSteadyStateAdaptationFromLegacyValidation(testCase)
%% Production-method version of isetvalidate/isetbio/cones/v_ibio_coneAdaptation.m.

outerSegment = osBioPhys('eccentricity', 0);
backgroundRates = logspace(0, 5, 50);
state = outerSegment.osAdaptSteadyState(backgroundRates);

testCase.verifySize(state.bgCur, size(backgroundRates));
testCase.verifyGreaterThan(state.bgCur, zeros(size(backgroundRates)));
testCase.verifyLessThanOrEqual(diff(state.bgCur), zeros(1, numel(backgroundRates)-1));
testCase.verifyEqual(size(state.opsin), size(backgroundRates));
testCase.verifyEqual(size(state.cGMP), size(backgroundRates));
end

function testSteadyStatePreservesInputShapeAndRepeatedValues(testCase)
outerSegment = osBioPhys('eccentricity', 15);
backgroundRates = [10 100; 10 100];
state = outerSegment.osAdaptSteadyState(backgroundRates);

testCase.verifySize(state.bgCur, size(backgroundRates));
testCase.verifyEqual(state.bgCur(1, :), state.bgCur(2, :), 'AbsTol', 1e-12);
testCase.verifyLessThan(state.bgCur(1, 2), state.bgCur(1, 1));
end

function testSteadyStateAdaptationGoldenValues(testCase)
%% Fingerprint foveal and peripheral steady-state current calculations.

currentTolerancePA = 1e-6;
backgroundRates = [1 10 100 1000 10000 100000];

fovealState = osBioPhys('eccentricity', 0).osAdaptSteadyState(backgroundRates);
peripheralState = osBioPhys('eccentricity', 15).osAdaptSteadyState(backgroundRates);

expectedFovealCurrent = [ ...
    86.1483675837284 86.1224309036816 85.8650666049445 ...
    83.4858806934175 69.5045417011094 38.2895494633141];
expectedPeripheralCurrent = [ ...
    86.1507790852277 86.1466704871920 86.1053863073252 ...
    85.6982063844960 82.0890342503994 64.3744300562738];

testCase.verifyEqual(fovealState.bgCur, expectedFovealCurrent, ...
    'AbsTol', currentTolerancePA);
testCase.verifyEqual(peripheralState.bgCur, expectedPeripheralCurrent, ...
    'AbsTol', currentTolerancePA);
end
