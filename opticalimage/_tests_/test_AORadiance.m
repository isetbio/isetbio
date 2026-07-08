function tests = test_AORadiance()
tests = functiontests(localfunctions);
end

function testOutputDimensions(testCase)
% Output must be a column vector the same length as wls, with exactly one
% nonzero element (at the component wavelength).
wls = (400:10:700)';
radiance = AOMonochromaticCornealPowerToRadiance(wls, 550, 1, 3, 1);

testCase.verifySize(radiance, [numel(wls) 1]);
testCase.verifyGreaterThan(radiance(wls == 550), 0);
testCase.verifyEqual(sum(radiance ~= 0), 1);
end

function testLinearInPower(testCase)
% Doubling corneal power must double radiance exactly (linear system).
wls = (400:10:700)';

r1 = AOMonochromaticCornealPowerToRadiance(wls, 550, 1, 3, 1);
r2 = AOMonochromaticCornealPowerToRadiance(wls, 550, 2, 3, 1);

testCase.verifyEqual(r2, 2 * r1, 'AbsTol', 1e-12);
end

function testInversePupilAreaScaling(testCase)
% Radiance is proportional to 1/pupilArea = 1/(pi*(d/2)^2), so
% r(2mm) / r(4mm) = (4/2)^2 = 4.
wls = (400:10:700)';

r2mm = AOMonochromaticCornealPowerToRadiance(wls, 550, 1, 2, 1);
r4mm = AOMonochromaticCornealPowerToRadiance(wls, 550, 1, 4, 1);

testCase.verifyEqual(r2mm, 4 * r4mm, 'RelTol', 1e-10);
end

function testWavelengthNotInWlsErrors(testCase)
% A component wavelength that does not appear in wls must raise an error.
wls = (400:10:700)';

localVerifyErrorContains(testCase, ...
    @() AOMonochromaticCornealPowerToRadiance(wls, 555, 1, 3, 1), ...
    'funky');
end

function testTwoComponentSuperposition(testCase)
% Radiance from two components equals the sum of each computed separately.
wls = (400:10:700)';

r500  = AOMonochromaticCornealPowerToRadiance(wls, 500, 1, 3, 1);
r600  = AOMonochromaticCornealPowerToRadiance(wls, 600, 1, 3, 1);
rBoth = AOMonochromaticCornealPowerToRadiance(wls, [500 600], [1 1], 3, 1);

testCase.verifyEqual(rBoth, r500 + r600, 'AbsTol', 1e-12);
end

% -------------------------------------------------------------------------
% Helper
% -------------------------------------------------------------------------

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
