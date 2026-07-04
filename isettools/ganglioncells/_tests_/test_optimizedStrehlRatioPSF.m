function tests = test_optimizedStrehlRatioPSF()
tests = functiontests(localfunctions);
end

function testAcceptsCanonicalZernikeDataBaseField(testCase)
opticsParams = localOpticsParams();
opticsParams.ZernikeDataBase = 'Polans2015';

localVerifyOptimizationUsesZernikeDataBase(testCase, opticsParams);
end

function testAcceptsLowercaseZernikeDataBaseField(testCase)
opticsParams = localOpticsParams();
opticsParams.zernikeDataBase = 'Polans2015';

localVerifyOptimizationUsesZernikeDataBase(testCase, opticsParams);
end

function localVerifyOptimizationUsesZernikeDataBase(testCase, opticsParams)
examinedDefocus = [-0.25 0 0.25];
fakeMosaic = FakeConeMosaicForStrehlTest('Polans2015');

[theOI, thePSF, optimalDefocus, optimalRatio, strehlRatio] = ...
    RGCMosaicConstructor.helper.optics.optimizedStrehlRatioPSF(...
    examinedDefocus, fakeMosaic, [0 0], opticsParams, 3, 1, false, false);

testCase.verifyEqual(optimalDefocus, 0.25);
testCase.verifyEqual(optimalRatio, 0.5, 'AbsTol', 1e-12);
testCase.verifyEqual(strehlRatio, [0.375 0.4375 0.5], 'AbsTol', 1e-12);
testCase.verifyEqual(theOI.refractiveErrorDiopters, optimalDefocus);
testCase.verifyEqual(thePSF.supportWavelength, [500 550 600]);
end

function opticsParams = localOpticsParams()
opticsParams = struct(...
    'subjectID', 1, ...
    'pupilDiameterMM', 3, ...
    'zeroCenterPSF', true, ...
    'noLCA', false);
end
