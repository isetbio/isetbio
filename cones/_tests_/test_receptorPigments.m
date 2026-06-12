function tests = test_receptorPigments()
tests = functiontests(localfunctions);
end

function testSharedSpectralContract(testCase)
wave = (500:10:540)';
absorbance = [0.1 0.2 0.3; 0.3 0.4 0.5; 1 1 1; 0.3 0.4 0.5; 0.1 0.2 0.3];
classes = {@cPhotoPigment, @photoPigment};

for ii = 1:numel(classes)
    pigment = classes{ii}('wave', wave, 'absorbance', absorbance);
    testCase.verifyEqual(pigment.wave, wave);
    testCase.verifySize(pigment.absorptance, size(absorbance));
    testCase.verifyGreaterThanOrEqual(pigment.quantalEfficiency, zeros(size(absorbance)));
    testCase.verifyLessThanOrEqual(pigment.quantalEfficiency, ones(size(absorbance)));
    testCase.verifyEqual(max(pigment.quantaFundamentals), ones(1, 3), 'AbsTol', 1e-12);
end
end

function testSubclassGeometry(testCase)
wave = [500 510];
absorbance = ones(2, 3);
diskPigment = cPhotoPigment('wave', wave, 'absorbance', absorbance, 'diameter', 4e-6);
rectPigment = photoPigment('wave', wave, 'absorbance', absorbance, ...
    'width', 4e-6, 'height', 5e-6, 'pdWidth', 2e-6, 'pdHeight', 3e-6);

testCase.verifyEqual(diskPigment.area, pi*(2e-6)^2, 'RelTol', 1e-12);
testCase.verifyEqual(rectPigment.area, 20e-12, 'RelTol', 1e-12);
testCase.verifyEqual(rectPigment.pdArea, 6e-12, 'RelTol', 1e-12);
testCase.verifyEqual([rectPigment.gapWidth rectPigment.gapHeight], [2e-6 2e-6], 'AbsTol', 1e-18);
end

function testDefaultReceptorPigmentGoldenValues(testCase)
%% Fingerprint the default Stockman-Sharpe receptor pigment data.

spectralTolerance = 1e-12;
pigment = cPhotoPigment('wave', [450 550 650]);

expectedAbsorbance = [ ...
    0.133838173147422 0.215010663808053 0.681239779965915; ...
    0.990923208090168 0.845162073638669 0.000262228565989288; ...
    0.0758280695937327 0.00599528730348394 3.29670434122775e-08];
expectedQuantalEfficiency = [ ...
    0.0952016434207283 0.146187429955574 0.310697133558703; ...
    0.453633542263389 0.414708926697902 0.000160994847231940; ...
    0.0557320888215983 0.00458570881485280 2.02425124070279e-08];

testCase.verifyEqual(pigment.absorbance, expectedAbsorbance, ...
    'AbsTol', spectralTolerance);
testCase.verifyEqual(pigment.quantalEfficiency, expectedQuantalEfficiency, ...
    'AbsTol', spectralTolerance);
end

function testWavelengthSetterResamplesSpectra(testCase)
classes = {@cPhotoPigment, @photoPigment};
newWave = [450 525 600 675]';

for ii = 1:numel(classes)
    pigment = classes{ii}();
    pigment.wave = newWave;

    testCase.verifyEqual(pigment.wave, newWave);
    testCase.verifySize(pigment.absorbance, [numel(newWave) 3]);
    testCase.verifySize(pigment.quantalEfficiency, [numel(newWave) 3]);
    testCase.verifyTrue(all(isfinite(pigment.quantalEfficiency(:))));
end
end

function testOpticalDensitySetterUpdatesAbsorptance(testCase)
wave = (500:10:540)';
absorbance = repmat([0.2 0.5 0.8], numel(wave), 1);
pigment = cPhotoPigment('wave', wave, ...
    'absorbance', absorbance, ...
    'opticalDensity', [0.2 0.2 0.2]);
lowDensityAbsorptance = pigment.absorptance;

pigment.opticalDensity = [0.4 0.4 0.4];

testCase.verifyGreaterThan(pigment.absorptance, lowDensityAbsorptance);
testCase.verifyEqual(pigment.absorptance, ...
    1-10.^(-absorbance*diag(pigment.opticalDensity)), ...
    'AbsTol', 1e-12);
end

function testRejectsInconsistentSpectralDimensions(testCase)
localVerifyErrorContains(testCase, ...
    @() cPhotoPigment('wave', [500 510], ...
    'absorbance', ones(2, 2), ...
    'opticalDensity', [0.5 0.5 0.4]), ...
    'optical density dimensionality does not match that of absorbance');

localVerifyErrorContains(testCase, ...
    @() cPhotoPigment('opticalDensity', [0.5 0.5 0.4], ...
    'peakEfficiency', [0.5 0.5]), ...
    'optical density dimensionality does not match that of peak efficiency');
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
