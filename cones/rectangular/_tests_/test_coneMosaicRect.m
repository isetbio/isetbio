function tests = test_coneMosaicRect()
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
scene = sceneCreate('uniform ee');
oi = oiCreate('human');
testCase.TestData.oi = oiCompute(oi, scene, 'pad value', 'mean');
end

function testSmallMosaicGeometry(testCase)
pattern = [2 3 4; 4 2 3];
mosaic = coneMosaicRect('pattern', pattern, 'noiseFlag', 'none', 'useParfor', false);

testCase.verifyEqual(mosaic.pattern, pattern);
testCase.verifyEqual(mosaic.mosaicSize, [2 3]);
testCase.verifyEqual([mosaic.rows mosaic.cols], [2 3]);
testCase.verifyEqual(mosaic.tSamples, 1);
testCase.verifySize(mosaic.coneLocs, [numel(pattern) 2]);
testCase.verifyEqual(sum(mosaic.spatialDensity), 1, 'AbsTol', 1e-12);
end

function testSettersMaintainDimensions(testCase)
mosaic = coneMosaicRect('pattern', ones(2, 3)*2, 'noiseFlag', 'none', 'useParfor', false);
mosaic.mosaicSize = [4 5];

testCase.verifyEqual(mosaic.mosaicSize, [4 5]);
testCase.verifySize(mosaic.pattern, [4 5]);

mosaic.absorptions = ones(4, 5, 3);
testCase.verifyEqual(mosaic.absorptions, ones(4, 5, 3));

mosaic.emPositions = zeros(3, 2);
testCase.verifyEqual(mosaic.tSamples, 3);
testCase.verifyEqual(mosaic.timeAxis, (0:2)*mosaic.integrationTime, 'AbsTol', 1e-12);
end

function testComputedAbsorptionDimensionsForSingleTimePoint(testCase)
mosaicSizes = {[1 1], [1 100], [20 20]};

for ii = 1:numel(mosaicSizes)
    mosaicSize = mosaicSizes{ii};
    mosaic = coneMosaicRect('size', mosaicSize, ...
        'noiseFlag', 'none', 'useParfor', false);

    absorptions = mosaic.compute(testCase.TestData.oi);

    testCase.verifyEqual(size(absorptions, 1), mosaicSize(1));
    testCase.verifyEqual(size(absorptions, 2), mosaicSize(2));
    testCase.verifyEqual(size(absorptions, 3), 1);
    testCase.verifyEqual(size(mosaic.absorptions, 1), mosaicSize(1));
    testCase.verifyEqual(size(mosaic.absorptions, 2), mosaicSize(2));
    testCase.verifyEqual(size(mosaic.absorptions, 3), 1);
end
end

function testComputedAbsorptionDimensionsForMultipleTimePoints(testCase)
testCases = {
    [1 1], 10
    [1 100], 10
    [20 20], 20
    };

for ii = 1:size(testCases, 1)
    mosaicSize = testCases{ii, 1};
    nTimePoints = testCases{ii, 2};
    mosaic = coneMosaicRect('size', mosaicSize, ...
        'noiseFlag', 'none', 'useParfor', false);
    mosaic.emGenSequence(nTimePoints, 'rSeed', 1, 'useParfor', false);

    absorptions = mosaic.compute(testCase.TestData.oi);

    testCase.verifySize(absorptions, [1 mosaicSize nTimePoints]);
    testCase.verifySize(mosaic.absorptions, [mosaicSize nTimePoints]);
    testCase.verifyEqual(mosaic.tSamples, nTimePoints);
end
end
