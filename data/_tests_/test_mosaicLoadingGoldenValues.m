function tests = test_mosaicLoadingGoldenValues()
% Golden-value regression tests for mosaic loading and source lattices.
%
% The values in this file are intentionally independent of the current MAT
% filenames. They protect the public loading contracts used by cone-mosaic
% examples and the standard source-lattice crops used to build cone and mRGC
% mosaics. Future SDR/cache loaders must preserve these results.

tests = functiontests(localfunctions);
end

function testConeMosaicExampleLoadingContracts(testCase)
% These are the two public forms used by the cone-mosaic examples.
byParameters = mosaicLoad([0.5 0.5], [0 0]);
byName = mosaicLoad('cmosaic_0.5-0.5_0.0-0.0.mat');

testCase.verifyEqual(byParameters.conesNum, 3822);
testCase.verifyEqual(byParameters.sizeDegs, ...
    [0.499915389531311 0.499838323158561], 'AbsTol', 1e-12);
testCase.verifyEqual(byParameters.eccentricityDegs, ...
    [-9.50412953337754e-06 2.04359107710966e-05], 'AbsTol', 1e-12);
testCase.verifyEqual(byParameters.coneRFpositionsDegs([1 end], :), ...
    [0.248109599164229 -0.134981796942376; ...
    -0.246855082615315 -0.12341241756038], 'AbsTol', 1e-12);

testCase.verifyEqual(byParameters.coneRFpositionsDegs, ...
    byName.coneRFpositionsDegs, 'AbsTol', 1e-12);
testCase.verifyEqual(byParameters.coneTypes, byName.coneTypes);
end

function testConeLatticeCropGoldenValues(testCase)
% Match the source lattice used by the default cMosaic construction path.
positionsMicrons = retinalattice.import.finalConePositions( ...
    58, [0 0], [0.5 0.5], 'right eye', []);

testCase.verifyEqual(size(positionsMicrons), [3525 2]);
testCase.verifyEqual(positionsMicrons([1 end], :), ...
    [66.8192901611328 55.5558128356934; ...
    -66.8488845825195 45.7090034484863], 'AbsTol', 1e-10);
end

function testMRGCLatticeCropGoldenValues(testCase)
% Match the 64-degree mRGC lattice used by standard circuit synthesis.
sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons( ...
    [0.5 0.5], 0);
[positionsMicrons, positionsDegs] = retinalattice.import.finalMRGCPositions( ...
    64, [0 0], sizeMicrons, 'right eye', ...
    @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x));

testCase.verifyEqual(size(positionsMicrons), [3265 2]);
testCase.verifyEqual(positionsMicrons([1 end], :), ...
    [66.8051986694336 39.4965705871582; ...
    -66.0864868164062 -64.4742279052734], 'AbsTol', 1e-10);
testCase.verifyEqual(positionsDegs([1 end], :), ...
    [0.237824562399033 0.1405428419368; ...
    -0.235263168875128 -0.229517512189418], 'AbsTol', 1e-12);
end
