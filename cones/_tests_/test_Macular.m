function tests = test_Macular()
tests = functiontests(localfunctions);
end

function testDependentSpectralProperties(testCase)
wave = (500:10:540)';
unitDensity = [0; 0.5; 1; 0.5; 0];
macular = Macular('wave', wave, 'unitDensity', unitDensity, 'density', 0.4);

testCase.verifyEqual(macular.wave, wave);
testCase.verifyEqual(macular.spectralDensity, 0.4*unitDensity, 'AbsTol', 1e-12);
testCase.verifyEqual(macular.transmittance, 10.^(-macular.spectralDensity), 'AbsTol', 1e-12);
testCase.verifyEqual(macular.absorptance + macular.transmittance, ones(size(wave)), 'AbsTol', 1e-12);
end

function testEccentricityDensity(testCase)
macular = Macular('wave', (500:10:540)', 'unitDensity', [0; 0.5; 1; 0.5; 0]);
ecc = [0 1 5];

testCase.verifyEqual(macular.eccDensity(ecc), macular.eccDensity([], 'eccDegs2', ecc.^2), 'AbsTol', 1e-12);
testCase.verifyEqual(macular.eccDensity(0), macular.density, 'AbsTol', 1e-12);
end

function testDensityScalingFromLegacyPigmentValidation(testCase)
%% Quantitative version of isetvalidate/isetbio/cones/v_ibio_pigments.m.

wave = (500:10:540)';
unitDensity = [0; 0.5; 1; 0.5; 0];
macular = Macular('wave', wave, 'unitDensity', unitDensity);
densityList = 0:0.1:0.5;

absorptance = zeros(numel(wave), numel(densityList));
for ii = 1:numel(densityList)
    macular.density = densityList(ii);
    testCase.verifyEqual(macular.spectralDensity, densityList(ii)*unitDensity, 'AbsTol', 1e-12);
    absorptance(:, ii) = macular.absorptance;
end

testCase.verifyGreaterThanOrEqual(diff(absorptance, 1, 2), zeros(numel(wave), numel(densityList)-1));
testCase.verifyEqual(absorptance(:, 1), zeros(numel(wave), 1), 'AbsTol', 1e-12);
end

function testDefaultMacularGoldenValues(testCase)
%% Fingerprint the default macular pigment data and density convention.

spectralTolerance = 1e-12;
macular = Macular('wave', [450 550 650]);

expectedUnitDensity = [0.946322067594433 0 0]';
expectedTransmittance = [0.466430859836272 1 1]';

testCase.verifyEqual(macular.unitDensity, expectedUnitDensity, ...
    'AbsTol', spectralTolerance);
testCase.verifyEqual(macular.transmittance, expectedTransmittance, ...
    'AbsTol', spectralTolerance);
end
