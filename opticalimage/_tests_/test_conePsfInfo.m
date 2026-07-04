function tests = test_conePsfInfo()
tests = functiontests(localfunctions);
end

function testDefaultStructureFields(testCase)
% Default construction must populate all four fields with consistent sizes.
info = conePsfInfoCreate;

testCase.verifyTrue(isfield(info, 'wavelengths'));
testCase.verifyTrue(isfield(info, 'spectralSensitivities'));
testCase.verifyTrue(isfield(info, 'spectralWeighting'));
testCase.verifyTrue(isfield(info, 'coneWeighting'));

nWave = numel(info.wavelengths);
testCase.verifySize(info.spectralSensitivities, [3, nWave]);
testCase.verifySize(info.spectralWeighting, [nWave, 1]);
testCase.verifySize(info.coneWeighting, [3, 1]);
end

function testDefaultWeightingsNormalize(testCase)
% Both default weight vectors must sum to exactly 1.
info = conePsfInfoCreate;

testCase.verifyEqual(sum(info.spectralWeighting), 1, 'AbsTol', 1e-12);
testCase.verifyEqual(sum(info.coneWeighting),     1, 'AbsTol', 1e-12);
end

function testGetNormalizesOnReturn(testCase)
% conePsfInfoGet divides by sum on every read, so a scaled field should
% still return a unit-sum vector.
info = conePsfInfoCreate;
info.spectralWeighting = 2 * info.spectralWeighting;  % bypass Set normalisation

retrieved = conePsfInfoGet(info, 'spectralWeighting');
testCase.verifyEqual(sum(retrieved), 1, 'AbsTol', 1e-12);
end

function testSetSpectralWeightingNormalizes(testCase)
% conePsfInfoSet normalises spectralWeighting on assignment.
info = conePsfInfoCreate;
nWave = numel(info.wavelengths);
unnormalized = 2 * ones(nWave, 1);  % sums to 2*nWave

info = conePsfInfoSet(info, 'spectralWeighting', unnormalized);
retrieved = conePsfInfoGet(info, 'spectralWeighting');

testCase.verifyEqual(sum(retrieved), 1, 'AbsTol', 1e-12);
end

function testSetWavelengthsResplinesSpectra(testCase)
% Changing wavelengths via Set must respline spectralSensitivities to the
% new sampling.
info    = conePsfInfoCreate;
newWave = (400:5:700)';  % 61 wavelengths (vs. default 31)

info = conePsfInfoSet(info, 'wavelengths', newWave);

testCase.verifyEqual(conePsfInfoGet(info, 'wavelengths'), newWave);
T = conePsfInfoGet(info, 'spectralSensitivities');
testCase.verifySize(T, [3, numel(newWave)]);
end

function testInconsistentSensitivityErrors(testCase)
% Setting spectralSensitivities with the wrong number of columns must error.
info = conePsfInfoCreate;
wrongCols = ones(3, 10);  % default has 31 wavelength columns

localVerifyErrorContains(testCase, ...
    @() conePsfInfoSet(info, 'spectralSensitivities', wrongCols), ...
    'not consistent');
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
