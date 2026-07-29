function tests = test_mRGCMosaicPrebakedGoldenValues()
% Golden-value regression tests for a prebaked mRGC mosaic.
%
% These run against the mosaic that ships in ganglioncells/mosaics, i.e. the
% artifact other scientists actually load, rather than against a stub. They pin
% both the structure of the mosaic and a deterministic noise-free response, so
% that a change anywhere in the connectivity or compute path shows up as a
% numeric failure instead of a plausible-looking but different answer.
%
% The mosaic used here is the smallest shipped one (0.5 x 0.5 deg at 2 deg
% nasal), chosen so the suite stays fast; it loads in well under a second.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.theMosaic = mRGCMosaic.loadPrebakedMosaic( ...
    'human', 'JCNpaperDefaultSubject', 'JCNpaperNasal2DegsTinyMosaic', ...
    'default', 'computeTheMosaicOptics', false);
end

% -------------------------------------------------------------------------
% Structure
% -------------------------------------------------------------------------

function testMosaicGeometry(testCase)
theMosaic = testCase.TestData.theMosaic;

testCase.verifyEqual(theMosaic.rgcsNum, 557, ...
    'Number of mRGCs in the prebaked mosaic changed.');
testCase.verifyEqual(theMosaic.eccentricityDegs, [2 0], 'AbsTol', 1e-12);
testCase.verifyEqual(theMosaic.sizeDegs, [0.5 0.5], 'AbsTol', 1e-12);
testCase.verifyEqual(theMosaic.whichEye, 'right eye');
end

function testInputConeMosaicComposition(testCase)
% Cone counts by type. The L:M ratio and S-cone fraction are the properties a
% user would notice first if the input mosaic were regenerated differently.
theConeMosaic = testCase.TestData.theMosaic.inputConeMosaic;

numL = sum(theConeMosaic.coneTypes == theConeMosaic.LCONE_ID);
numM = sum(theConeMosaic.coneTypes == theConeMosaic.MCONE_ID);
numS = sum(theConeMosaic.coneTypes == theConeMosaic.SCONE_ID);

testCase.verifyEqual(numel(theConeMosaic.coneTypes), 3765, ...
    'Number of input cones changed.');
testCase.verifyEqual([numL numM numS], [2270 1118 377], ...
    'Cone type counts changed.');

% Sanity bounds independent of the exact golden numbers.
testCase.verifyGreaterThan(numL/numM, 1.5, 'Expected an L-dominant mosaic.');
testCase.verifyLessThan(numL/numM, 2.5, 'L:M ratio outside the expected range.');
testCase.verifyLessThan(numS/(numL+numM+numS), 0.15, ...
    'S-cone fraction is higher than expected for human retina.');
end

function testRFcenterConnectivity(testCase)
% At 2 deg eccentricity midget RGC centers should draw from one or two cones.
theMosaic = testCase.TestData.theMosaic;
centerMatrix = theMosaic.rgcRFcenterConeConnectivityMatrix;

testCase.verifyEqual(size(centerMatrix), [3765 557]);
testCase.verifyEqual(nnz(centerMatrix), 697, ...
    'Number of cone-to-RF-center connections changed.');
testCase.verifyEqual(full(sum(centerMatrix(:))), 697, 'AbsTol', 1e-12, ...
    'RF center weights are no longer unit-valued.');

conesPerCenter = full(sum(centerMatrix > 0, 1));
testCase.verifyEqual(min(conesPerCenter), 1);
testCase.verifyEqual(max(conesPerCenter), 2);
testCase.verifyEqual(mean(conesPerCenter), 1.251346499102, 'RelTol', 1e-10, ...
    'Mean cones per RF center changed.');

testCase.verifyFalse(any(isnan(nonzeros(centerMatrix))), ...
    'RF center connectivity must not contain NaN weights.');
end

% -------------------------------------------------------------------------
% Surround integrity
% -------------------------------------------------------------------------

function testSurroundConnectivityKnownDegenerateCells(testCase)
% KNOWN DEFECT, pinned deliberately.
%
% Four of the 557 RGCs in this mosaic carry all-NaN surround weights, so their
% surround sums to zero and they behave as center-only cells: they respond at
% the full center amplitude while their neighbours respond at roughly a tenth
% of it. Across all 16 shipped mosaics this affects 15 of 102619 RGCs.
%
% This test documents the defect so it cannot get quietly worse. When the
% mosaics are regenerated without the degenerate cells, this test should fail
% and be tightened to require zero such cells.
theMosaic = testCase.TestData.theMosaic;
surroundMatrix = theMosaic.rgcRFsurroundConeConnectivityMatrix;

degenerateRGCindices = find(any(isnan(surroundMatrix), 1));
testCase.verifyEqual(degenerateRGCindices, [156 264 282 303], ...
    'The set of RGCs with NaN surround weights changed.');

% Every other RGC must have a real, positive surround.
surroundWeightSums = full(sum(surroundMatrix, 1, 'omitnan'));
healthyRGCindices = setdiff(1:theMosaic.rgcsNum, degenerateRGCindices);
testCase.verifyTrue(all(surroundWeightSums(healthyRGCindices) > 0), ...
    'Every non-degenerate RGC must pool a non-zero surround.');
testCase.verifyEqual(surroundWeightSums(degenerateRGCindices), zeros(1,4), ...
    'AbsTol', 1e-12);
end

% -------------------------------------------------------------------------
% Deterministic response
% -------------------------------------------------------------------------

function testNoiseFreeResponseToSpatialModulation(testCase)
% End-to-end golden values for the noise-free compute path.
%
% The drive is an analytic function of cone position rather than a random
% draw, so the expected values do not depend on the RNG implementation.
theMosaic = testCase.TestData.theMosaic;

coneModulation = localConeModulation(theMosaic);
testCase.verifyEqual(sum(coneModulation), -0.885408209986, 'RelTol', 1e-9, ...
    'The test drive itself changed; expected responses below are stale.');

theMosaic.noiseFlag = 'none';
theResponse = squeeze(theMosaic.compute( ...
    reshape(coneModulation, [1 1 numel(coneModulation)]), 0));

testCase.verifyEqual(numel(theResponse), 557);
testCase.verifyFalse(any(isnan(theResponse)), ...
    'Noise-free responses must not contain NaN.');

% Evaluate the healthy cells; the degenerate four are covered above.
degenerateRGCindices = [156 264 282 303];
healthyResponse = theResponse(setdiff(1:theMosaic.rgcsNum, degenerateRGCindices));

testCase.verifyEqual(mean(healthyResponse), 0.000034610992, 'AbsTol', 1e-9, ...
    'Mean noise-free response changed.');
testCase.verifyEqual(std(healthyResponse), 0.033861169695, 'RelTol', 1e-8, ...
    'Response standard deviation changed.');
testCase.verifyEqual(min(healthyResponse), -0.085553085396, 'RelTol', 1e-8);
testCase.verifyEqual(max(healthyResponse), 0.085732251675, 'RelTol', 1e-8);

% A center-surround mosaic driven by a zero-mean grating should produce a
% near-zero mean response; this is the qualitative check behind the numbers.
testCase.verifyLessThan(abs(mean(healthyResponse)), 0.01 * std(healthyResponse), ...
    'Mean response should be small relative to its spread.');
end

function testDegenerateCellsRespondWithoutSurround(testCase)
% Companion to testSurroundConnectivityKnownDegenerateCells: shows the defect
% is not merely cosmetic. Under a spatially uniform drive a healthy
% center-surround cell largely cancels, while the four degenerate cells pass
% the drive through at full amplitude.
theMosaic = testCase.TestData.theMosaic;
theMosaic.noiseFlag = 'none';

uniformDrive = 0.1;
numCones = size(theMosaic.rgcRFcenterConeConnectivityMatrix, 1);
theResponse = squeeze(theMosaic.compute( ...
    repmat(uniformDrive, [1 1 numCones]), 0));

degenerateRGCindices = [156 264 282 303];
healthyRGCindices = setdiff(1:theMosaic.rgcsNum, degenerateRGCindices);

testCase.verifyEqual(theResponse(degenerateRGCindices), ...
    repmat(uniformDrive, [4 1]), 'AbsTol', 1e-12, ...
    'Degenerate cells should pass the uniform drive through unattenuated.');

testCase.verifyLessThan(max(theResponse(healthyRGCindices)), 0.5 * uniformDrive, ...
    'Healthy cells should attenuate a uniform drive via their surround.');
end

% -------------------------------------------------------------------------

function coneModulation = localConeModulation(theMosaic)
% Deterministic, RNG-free drive: a 0.1 deg period grating across the mosaic.
conePositionsDegs = theMosaic.inputConeMosaic.coneRFpositionsDegs;
xRelativeDegs = conePositionsDegs(:,1) - theMosaic.eccentricityDegs(1);
coneModulation = 0.1 * sin(2*pi*xRelativeDegs/0.1);
end
