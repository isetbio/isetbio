function tests = test_cMosaicCompute()
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
scene = sceneCreate('uniform ee');
scene = sceneSet(scene, 'fov', 0.5);
oi = oiCreate('human');
testCase.TestData.oi = oiCompute(oi, scene, 'pad value', 'mean');
end

function testDeterministicComputeContract(testCase)
mosaic = localAllLMosaic(5/1000);

responses = mosaic.compute(testCase.TestData.oi, ...
    'nTrials', 2, 'nTimePoints', 3);

% Without noise or eye movements, compute returns one response trial even
% when more trials are requested.
testCase.verifySize(responses, [1 3 mosaic.conesNum]);
testCase.verifyEqual(responses(1, 1, :), responses(1, 3, :), ...
    'AbsTol', 1e-12);
testCase.verifyGreaterThanOrEqual(responses, zeros(size(responses)));
testCase.verifyTrue(all(isfinite(responses(:))));
end

function testUniformResponseAndIntegrationTimeScaling(testCase)
shortMosaic = localAllLMosaic(5/1000);
longMosaic = localAllLMosaic(10/1000);

shortResponse = squeeze(shortMosaic.compute(testCase.TestData.oi));
longResponse = squeeze(longMosaic.compute(testCase.TestData.oi));

testCase.verifyLessThan(range(shortResponse)/mean(shortResponse), 2e-5);
testCase.verifyEqual(longResponse, 2*shortResponse, 'RelTol', 1e-10);
end

function testSeededEyeMovementsAndResponseDimensions(testCase)
firstMosaic = localAllLMosaic(5/1000);
secondMosaic = localAllLMosaic(5/1000);

durationSeconds = 20/1000;
firstMosaic.emGenSequence(durationSeconds, ...
    'nTrials', 2, 'randomSeed', 3, 'useParfor', false);
secondMosaic.emGenSequence(durationSeconds, ...
    'nTrials', 2, 'randomSeed', 3, 'useParfor', false);

testCase.verifyEqual(firstMosaic.fixEMobj.emPosMicrons, ...
    secondMosaic.fixEMobj.emPosMicrons, 'AbsTol', 1e-12);

responses = firstMosaic.compute(testCase.TestData.oi, ...
    'withFixationalEyeMovements', true);
testCase.verifySize(responses, [2 4 firstMosaic.conesNum]);
testCase.verifyTrue(all(isfinite(responses(:))));
end

function testComputeRequiresGeneratedEyeMovements(testCase)
mosaic = localAllLMosaic(5/1000);

localVerifyErrorContains(testCase, ...
    @() mosaic.compute(testCase.TestData.oi, ...
    'withFixationalEyeMovements', true), ...
    'cMosaic.emGenSequence() has not been called');
end

function testImportedDegreeAndMicronGeometryAgree(testCase)
degreeData = localConeData('degrees');
micronData = degreeData;
micronData.positionUnits = 'microns';
micronData.positions = 300*degreeData.positions;
micronData.lightGatheringApertureDiameters = ...
    300*degreeData.lightGatheringApertureDiameters;

degreeMosaic = localImportedMosaic(degreeData, 5/1000);
micronMosaic = localImportedMosaic(micronData, 5/1000);

testCase.verifyEqual(degreeMosaic.coneRFpositionsDegs, ...
    micronMosaic.coneRFpositionsDegs, 'AbsTol', 1e-12);
testCase.verifyEqual(degreeMosaic.coneRFpositionsMicrons, ...
    micronMosaic.coneRFpositionsMicrons, 'AbsTol', 1e-12);
testCase.verifyEqual(degreeMosaic.coneApertureDiametersDegs, ...
    micronMosaic.coneApertureDiametersDegs, 'AbsTol', 1e-12);
end

function testMalformedImportedConeData(testCase)
coneData = localConeData('degrees');

missingTypes = rmfield(coneData, 'types');
localVerifyErrorContains(testCase, ...
    @() localImportedMosaic(missingTypes, 5/1000), ...
    'expected a ''type'' field');

badPositions = coneData;
badPositions.positions = [badPositions.positions zeros(9, 1)];
localVerifyErrorContains(testCase, ...
    @() localImportedMosaic(badPositions, 5/1000), ...
    'positions'' should be an N x 2 matrix');

badApertures = coneData;
badApertures.lightGatheringApertureDiameters = ones(1, 8);
localVerifyErrorContains(testCase, ...
    @() localImportedMosaic(badApertures, 5/1000), ...
    'dimensions of ''positions''');
end

function testReassignConeTypesMaintainsPartition(testCase)
mosaic = localImportedMosaic(localConeData('degrees'), 5/1000);

mosaic.reassignTypeOfCones([], cMosaic.LCONE_ID);

testCase.verifyEqual(mosaic.lConeIndices, (1:mosaic.conesNum)');
testCase.verifyEmpty(mosaic.mConeIndices);
testCase.verifyEmpty(mosaic.sConeIndices);
testCase.verifyEmpty(mosaic.kConeIndices);
testCase.verifyEqual(mosaic.coneTypes, ...
    repmat(cMosaic.LCONE_ID, mosaic.conesNum, 1));
testCase.verifyEqual(mosaic.coneDensities, [1 0 0 0], 'AbsTol', 1e-12);
end

function mosaic = localAllLMosaic(integrationTime)
coneData = localConeData('degrees');
coneData.types = repmat(cMosaic.LCONE_ID, 1, size(coneData.positions, 1));
mosaic = localImportedMosaic(coneData, integrationTime);
end

function mosaic = localImportedMosaic(coneData, integrationTime)
mosaic = cMosaic( ...
    'coneData', coneData, ...
    'micronsPerDegree', 300, ...
    'integrationTime', integrationTime, ...
    'noiseFlag', 'none', ...
    'randomSeed', 1, ...
    'eccVaryingConeAperture', false, ...
    'eccVaryingOuterSegmentLength', false, ...
    'eccVaryingMacularPigmentDensity', false, ...
    'eccVaryingMacularPigmentDensityDynamic', false, ...
    'eccVaryingConeBlur', false, ...
    'useParfor', false);
end

function coneData = localConeData(positionUnits)
positions = [
    -0.1 -0.1
     0.0 -0.1
     0.1 -0.1
    -0.1  0.0
     0.0  0.0
     0.1  0.0
    -0.1  0.1
     0.0  0.1
     0.1  0.1
    ];
coneData = struct( ...
    'positionUnits', positionUnits, ...
    'positions', positions, ...
    'types', repmat([cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID], 1, 3), ...
    'lightGatheringApertureDiameters', 0.03*ones(1, size(positions, 1)), ...
    'blurApertureDiameterMicrons', 9);
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
