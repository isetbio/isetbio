function tests = test_cMosaicIntegrationFullOnly()
tests = functiontests(localfunctions);
end

function setup(testCase)
testCase.TestData.existingFigures = findall(groot, 'Type', 'figure');
end

function teardown(testCase)
allFigures = findall(groot, 'Type', 'figure');
testFigures = setdiff(allFigures, testCase.TestData.existingFigures);
testFigures = testFigures(ishghandle(testFigures));
if ~isempty(testFigures), close(testFigures); end
end

function testBundledLatticeAllLConstructorAndHarmonicResponse(testCase)
mosaic = cMosaic( ...
    'sizeDegs', [0.05 0.05], ...
    'eccentricityDegs', [0 0], ...
    'coneDensities', [1 0 0 0], ...
    'integrationTime', 10/1000, ...
    'noiseFlag', 'none', ...
    'randomSeed', 1, ...
    'useParfor', false);

testCase.verifyEqual(mosaic.lConeIndices, (1:mosaic.conesNum)');
testCase.verifyEmpty(mosaic.mConeIndices);
testCase.verifyEmpty(mosaic.sConeIndices);
testCase.verifyEmpty(mosaic.kConeIndices);
testCase.verifyEqual(mosaic.coneDensities, [1 0 0 0], 'AbsTol', 1e-12);
testCase.verifyTrue(all(isfinite(mosaic.coneRFpositionsDegs(:))));

responses = zeros(mosaic.conesNum, 2);
phases = [0 pi/2];
for ii = 1:numel(phases)
    params = harmonicP('freq', 4, 'contrast', 1, 'ph', phases(ii), ...
        'ang', 0, 'row', 128, 'col', 128);
    scene = sceneSet(sceneCreate('harmonic', params), 'fov', 0.08);
    oi = oiCompute(oiCreate('human'), scene, 'pad value', 'mean');
    responses(:, ii) = squeeze(mosaic.compute(oi));
end

testCase.verifyTrue(all(isfinite(responses(:))));
testCase.verifyGreaterThan(norm(responses(:, 1)-responses(:, 2)), 0);
end

function testDeadLeavesExampleSmoke(testCase)
result = localRunDeadLeavesExample();

testCase.verifyClass(result.cm, 'cMosaic');
testCase.verifySize(result.excitations, ...
    [1 numel(result.timeAxis) result.cm.conesNum]);
end

function testHarmonicExampleSmoke(testCase)
result = localRunHarmonicExample();

testCase.verifyClass(result.cm, 'cMosaic');
testCase.verifyNumElements(result.selectedConeIndices, 5);
testCase.verifySize(result.selectedResponses, [numel(result.timeAxis) 5]);
end

function testQuadratureExampleSmoke(testCase)
result = localRunQuadratureExample();

testCase.verifyClass(result.cm, 'cMosaic');
testCase.verifyEqual(result.cm.coneDensities, [1 0 0 0], 'AbsTol', 1e-12);
testCase.verifyTrue(all(isfinite(result.coneQuadratureEnergy)));
testCase.verifyLessThan( ...
    range(result.quadratureEnergy)/mean(result.quadratureEnergy), 1e-12);
end

function result = localRunDeadLeavesExample()
run(fullfile(isetbioRootPath, 'examples', 'conemosaic', 's_cmDeadLeaves.m'));
result = struct('cm', cm, 'excitations', excitations, 'timeAxis', timeAxis);
end

function result = localRunHarmonicExample()
run(fullfile(isetbioRootPath, 'examples', 'conemosaic', 's_cmHarmonic.m'));
result = struct( ...
    'cm', cm, ...
    'selectedConeIndices', selectedConeIndices, ...
    'selectedResponses', selectedResponses, ...
    'timeAxis', timeAxis);
end

function result = localRunQuadratureExample()
run(fullfile(isetbioRootPath, 'examples', 'conemosaic', 's_cmQuadrature.m'));
result = struct( ...
    'cm', cm, ...
    'coneQuadratureEnergy', coneQuadratureEnergy, ...
    'quadratureEnergy', quadratureEnergy);
end
