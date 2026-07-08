function tests = test_Lens()
tests = functiontests(localfunctions);
end

function testDependentSpectralProperties(testCase)
% Verify the three computed properties satisfy their defining formulas.
wave = (500:10:540)';
unitDensity = [0.2; 0.5; 1.0; 0.5; 0.2];
lens = Lens('wave', wave, 'unitDensity', unitDensity, 'density', 1.5);

testCase.verifyEqual(lens.wave, wave);
testCase.verifyEqual(lens.spectralDensity, 1.5 * unitDensity, 'AbsTol', 1e-12);
testCase.verifyEqual(lens.transmittance, 10.^(-lens.spectralDensity), 'AbsTol', 1e-12);
testCase.verifyEqual(lens.absorptance + lens.transmittance, ones(size(wave)), 'AbsTol', 1e-12);
end

function testDensitySetterUpdatesTransmittance(testCase)
% Changing density should update spectralDensity and transmittance consistently.
wave = (500:10:540)';
unitDensity = [0.2; 0.5; 1.0; 0.5; 0.2];
lens = Lens('wave', wave, 'unitDensity', unitDensity);

for d = [0.5 1.0 1.5 2.0]
    lens.density = d;
    testCase.verifyEqual(lens.spectralDensity, d * unitDensity, 'AbsTol', 1e-12);
    testCase.verifyEqual(lens.transmittance, 10.^(-lens.spectralDensity), 'AbsTol', 1e-12);
end
end

function testWavelengthSetterResamplesUnitDensity(testCase)
% Changing wave should resample unitDensity_ and keep the photon budget balanced.
lens = Lens();
newWave = [450 500 550 600 650]';

lens.wave = newWave;

testCase.verifyEqual(lens.wave, newWave);
testCase.verifySize(lens.unitDensity, size(newWave));
testCase.verifyTrue(all(isfinite(lens.transmittance)));
testCase.verifyEqual(lens.absorptance + lens.transmittance, ones(size(newWave)), 'AbsTol', 1e-12);
end

function testDefaultLensGoldenValues(testCase)
% Reproducibility: two default instances must agree to floating-point precision.
% Structural: default data must obey physical bounds and the known UV-heavy
% absorption of the human lens.
lens1 = Lens('wave', [450 550 650]);
lens2 = Lens('wave', [450 550 650]);

testCase.verifyEqual(lens1.unitDensity,   lens2.unitDensity,   'AbsTol', 1e-14);
testCase.verifyEqual(lens1.transmittance, lens2.transmittance, 'AbsTol', 1e-14);

% Human lens absorbs more at short wavelengths.
testCase.verifyGreaterThan(lens1.unitDensity(1), lens1.unitDensity(3));
testCase.verifyGreaterThan(lens1.unitDensity(1), 0);
testCase.verifyTrue(all(lens1.unitDensity   >= 0));
testCase.verifyTrue(all(lens1.transmittance >  0));
testCase.verifyTrue(all(lens1.transmittance <= 1));
end

function testLegacyGetWrapperConsistency(testCase)
% get(lens, param) must return the same values as direct property access.
lens = Lens('wave', (450:10:650)', 'density', 1.2);

testCase.verifyEqual(get(lens, 'absorbance'),    lens.unitDensity,    'AbsTol', 1e-15);
testCase.verifyEqual(get(lens, 'transmittance'), lens.transmittance,  'AbsTol', 1e-15);
testCase.verifyEqual(get(lens, 'absorptance'),   lens.absorptance,    'AbsTol', 1e-15);
testCase.verifyEqual(get(lens, 'spectraldensity'), lens.spectralDensity, 'AbsTol', 1e-15);
testCase.verifyEqual(get(lens, 'density'),       lens.density,        'AbsTol', 1e-15);
testCase.verifyEqual(get(lens, 'wave'),          lens.wave);
end
