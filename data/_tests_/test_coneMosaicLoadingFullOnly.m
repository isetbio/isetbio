function tests = test_coneMosaicLoadingFullOnly()
% Live SDR regression test for the public pre-baked cone mosaic loader.
%
% The FullOnly name excludes this test from the ordinary offline unit suite.

tests = functiontests(localfunctions);
end

function testConeMosaicExampleLoadingContracts(testCase)
% These are the two public forms used by cone-mosaic examples.
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
