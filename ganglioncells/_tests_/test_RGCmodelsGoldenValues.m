function tests = test_RGCmodelsGoldenValues()
% Golden-value regression tests for the analytic RGC models.
%
% These pin the published formulae that the mRGC mosaics are built on, so a
% change to a coefficient shows up as a test failure rather than as a silent
% shift in every synthesized mosaic. Expected values are quoted to more digits
% than the tolerance requires so the intended number is unambiguous.

tests = functiontests(localfunctions);
end

% -------------------------------------------------------------------------
% Watson (2014) retinal distance conversions, Appendix Equations A5 and A6.
% -------------------------------------------------------------------------

function testWatsonRhoDegsToMMs(testCase)
% Equation A5. The linear coefficient is the well known ~0.268 mm/deg near
% the fovea, so the 1 deg value doubles as a sanity check on the formula.
eccDegs = [1 5 10];
expectedMMs = [0.2683343691 1.3475261375 2.7059391000];

testCase.verifyEqual( ...
    RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs), expectedMMs, ...
    'RelTol', 1e-10, 'Watson Eq A5 (degs -> mm) changed.');
end

function testWatsonRhoMMsToDegs(testCase)
% Equation A6.
eccMMs = [0.2683343691 1.3475261375 2.7059391000];
expectedDegs = [0.958371581956 4.883619410680 9.931576670812];

testCase.verifyEqual( ...
    RGCmodels.Watson.convert.rhoMMsToDegs(eccMMs), expectedDegs, ...
    'RelTol', 1e-8, 'Watson Eq A6 (mm -> degs) changed.');
end

function testWatsonConversionsAreNotExactInverses(testCase)
% Watson's A5 and A6 are separately fitted polynomials, not an analytic
% inverse pair, so a degs -> mm -> degs round trip does not return the input.
% The discrepancy is largest in the parafovea, which is where most of the
% prebaked mRGC mosaics sit, so it is pinned here deliberately: code that
% needs an invertible mapping must not assume these two functions provide it.
%
% If a future change makes them true inverses, this test should fail and be
% updated rather than the behavior changing unnoticed.
eccDegs = [0.5 1 2 5 10];
roundTripped = RGCmodels.Watson.convert.rhoMMsToDegs( ...
    RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs));

relativeError = (roundTripped - eccDegs) ./ eccDegs;
expectedRelativeError = [-0.0442706 -0.0416284 -0.0365802 -0.0232766 -0.0068420];

testCase.verifyEqual(relativeError, expectedRelativeError, 'AbsTol', 1e-6, ...
    'Watson degs->mm->degs round-trip error changed.');

% Guard the qualitative fact rather than only the numbers.
testCase.verifyLessThan(abs(relativeError), 0.05, ...
    'Round-trip error should stay within 5% over the parafoveal range.');
end

% -------------------------------------------------------------------------
% Croner & Kaplan (1995) center/surround relationships.
% -------------------------------------------------------------------------

function testCronerKaplanSurroundToCenterRcRatio(testCase)
% Croner & Kaplan report a surround characteristic radius about 6.7 times the
% center radius. This single constant sets the surround extent of every
% synthesized mosaic.
testCase.verifyEqual( ...
    RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio, 6.7, ...
    'AbsTol', 1e-12, 'Croner & Kaplan surround/center Rc ratio changed.');
end

function testCronerKaplanSurroundCharacteristicRadius(testCase)
% Surround Rc in degrees, from the fit to the pooled P and M cell data.
% Returned as a column vector.
eccDegs = [1 5 10];
expectedRcDegs = [0.203000000000; 0.433920211926; 0.601860772073];

testCase.verifyEqual( ...
    RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(eccDegs), ...
    expectedRcDegs, 'RelTol', 1e-8, ...
    'Croner & Kaplan surround characteristic radius changed.');

% Surround radius must grow with eccentricity.
testCase.verifyTrue(all(diff(expectedRcDegs) > 0), ...
    'Surround radius should increase with eccentricity.');
end

function testCronerKaplanIntegratedSensitivityRatio(testCase)
% Surround/center integrated sensitivity for P cells. Croner & Kaplan report
% values near 0.5, i.e. the surround cancels roughly half the center.
eccDegs = [1 5 10];
expectedRatio = [0.473 0.501 0.536];

testCase.verifyEqual( ...
    RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(eccDegs), ...
    expectedRatio, 'RelTol', 1e-8, ...
    'Croner & Kaplan integrated sensitivity ratio changed.');

testCase.verifyTrue(all(expectedRatio > 0 & expectedRatio < 1), ...
    'Integrated sensitivity ratio must stay between 0 and 1.');
end
