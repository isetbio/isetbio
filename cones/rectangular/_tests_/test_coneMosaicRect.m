function tests = test_coneMosaicRect()
tests = functiontests(localfunctions);
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
