function tests = test_coneRectWindowFullOnly()
tests = functiontests(localfunctions);
end

function testDefaultPlotTypeWithData(testCase)
mosaic = coneMosaicRect('pattern', [2 3; 4 2], ...
    'noiseFlag', 'none', 'useParfor', false);
mosaic.absorptions = ones(2, 2, 2, 'single');
mosaic.current = -ones(2, 2, 2, 'single');

app = coneRectWindow(mosaic);
cleanup = onCleanup(@() delete(app)); %#ok<NASGU>

testCase.verifyClass(app, 'coneRectWindow_App');
testCase.verifyEqual(app.popupImageType.Value, 'Mean absorptions');
end

function testFailedComputePreservesAppAndResponses(testCase)
previousInitClear = ieSessionGet('init clear');
restoreInitClear = onCleanup( ...
    @() ieSessionSet('init clear', previousInitClear)); %#ok<NASGU>
ieSessionSet('init clear', false);
ieInit;

badOI = oiCreate('empty');
badOI.data.photons = ones(2, 2, 1);
ieAddObject(badOI);

mosaic = coneMosaicRect('pattern', [2 3; 4 2], ...
    'noiseFlag', 'none', 'useParfor', false);
previousAbsorptions = ones(2, 2, 2, 'single');
previousCurrent = -ones(2, 2, 2, 'single');
mosaic.absorptions = previousAbsorptions;
mosaic.current = previousCurrent;
previousAbsorptions = mosaic.absorptions;
previousCurrent = mosaic.current;

app = coneRectWindow(mosaic);
cleanup = onCleanup(@() delete(app)); %#ok<NASGU>

lastwarn('');
app.btnComputeImage.ButtonPushedFcn(app.btnComputeImage, []);
[~, warningID] = lastwarn;

testCase.verifyTrue(isvalid(app));
testCase.verifyEqual(mosaic.absorptions, previousAbsorptions);
testCase.verifyEqual(mosaic.current, previousCurrent);
testCase.verifyEqual(warningID, 'coneRectWindow:ComputeAbsorptionsFailed');
end
