function tests = test_cMosaic()
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
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
coneData = struct(...
    'positionUnits', 'degrees', ...
    'positions', positions, ...
    'types', repmat([cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID], 1, 3), ...
    'lightGatheringApertureDiameters', 0.03*ones(1, size(positions, 1)));
testCase.TestData.mosaic = cMosaic(...
    'coneData', coneData, ...
    'micronsPerDegree', 300, ...
    'noiseFlag', 'none', ...
    'randomSeed', 1, ...
    'useParfor', false);
testCase.TestData.excitations = reshape(1:54, [2 3 9]);
end

function testOuterSegmentLengthModel(testCase)
[lengths, fovealLength] = cMosaic.outerSegmentLengthFromEccentricity([0 4.86 39.90]);

testCase.verifyEqual(fovealLength, 47.81, 'AbsTol', 1e-12);
testCase.verifyEqual(lengths, [47.81 21.2 13.22], 'AbsTol', 1e-12);
testCase.verifyGreaterThanOrEqual(lengths(1:end-1), lengths(2:end));
end

function testNoisyInstancesFrozenSeed(testCase)
meanAbsorptions = reshape([1 10 30 100], [2 2]);
first = cMosaic.noisyInstances(meanAbsorptions, 'noiseFlag', 'frozen', 'seed', 9);
second = cMosaic.noisyInstances(meanAbsorptions, 'noiseFlag', 'frozen', 'seed', 9);

testCase.verifyEqual(first, second);
testCase.verifyGreaterThanOrEqual(first, zeros(size(first)));
end

function testPlotSelectsTrialAndTimePoint(testCase)
mosaic = testCase.TestData.mosaic;
excitations = testCase.TestData.excitations;

[plotData, figureHandle] = mosaic.plot('excitations', excitations, ...
    'trial', 2, 'time point', 'last', 'data only', true);

testCase.verifyEmpty(figureHandle);
testCase.verifyEqual(plotData.excitations, squeeze(excitations(2, end, :)));
testCase.verifyEqual(plotData.trial, 2);
testCase.verifyEqual(plotData.timePoint, size(excitations, 2));
end

function testPlotUsesSuppliedAxes(testCase)
mosaic = testCase.TestData.mosaic;
figureHandle = figure('Visible', 'off');
axesHandle = axes('Parent', figureHandle);

[plotData, returnedFigure] = mosaic.plot('excitations', ...
    testCase.TestData.excitations, 'axes handle', axesHandle);

testCase.verifyEqual(returnedFigure, figureHandle);
testCase.verifyEqual(plotData.figureHandle, figureHandle);
testCase.verifyEqual(plotData.axesHandle, axesHandle);
testCase.verifyNotEmpty(axesHandle.Children);
end

function testPlotProfileReturnsTwoDimensionalPositions(testCase)
mosaic = testCase.TestData.mosaic;

[plotData, figureHandle] = mosaic.plot('excitations horizontal line', ...
    testCase.TestData.excitations, 'time point', 1, ...
    'thickness', 0.04, 'data only', true);

testCase.verifyEmpty(figureHandle);
testCase.verifyNumElements(plotData.positions, 3);
for ii = 1:numel(plotData.positions)
    testCase.verifyEqual(size(plotData.positions{ii}, 2), 2);
end
end

function testPlotRequiresROI(testCase)
mosaic = testCase.TestData.mosaic;

testCase.verifyError(@() mosaic.plot('roi', testCase.TestData.excitations), ...
    'cMosaic:MissingROI');
end

function testPlotEyeMovements(testCase)
mosaic = testCase.TestData.mosaic;
mosaic.emGenSequence(20/1000, 'nTrials', 1, 'randomSeed', 1, 'useParfor', false);

[plotData, figureHandle] = mosaic.plot('excitations and eye movements', ...
    testCase.TestData.excitations, 'time point', 'last', 'data only', true);

testCase.verifyEmpty(figureHandle);
testCase.verifyEqual(plotData.timePoints, 1:size(testCase.TestData.excitations, 2));
testCase.verifySize(plotData.currentEMposition, [2 1]);
end
