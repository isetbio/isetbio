function tests = test_oisCreateFullOnly()
% Smoke tests for oisCreate. Named FullOnly because all paths call
% oiCreate('wvf human') and oiCompute, which are slow.
tests = functiontests(localfunctions);
end

function testHarmonicOisCreate(testCase)
hparams(2) = harmonicP;
hparams(2).freq = 6;
hparams(2).GaborFlag = 0.2;
hparams(1) = hparams(2);
hparams(1).contrast = 0;
sparams.fov = 0.5;  % small fov keeps compute fast
stimWeights = ones(1, 5) / 5;

ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);

testCase.verifyClass(ois, 'oiSequence');
testCase.verifyEqual(ois.length, numel(stimWeights));
frame = ois.frameAtIndex(1);
photons = oiGet(frame, 'photons');
testCase.verifyTrue(all(isfinite(photons(:))));
end

function testImpulseOisCreate(testCase)
sparams.fov = 0.5;
sparams.meanluminance = 100;
stimWeights = zeros(1, 5);
stimWeights(2) = 1;

ois = oisCreate('impulse', 'add', stimWeights, 'sceneParameters', sparams);

testCase.verifyClass(ois, 'oiSequence');
testCase.verifyEqual(ois.length, numel(stimWeights));
frame = ois.frameAtIndex(1);
photons = oiGet(frame, 'photons');
testCase.verifyTrue(all(isfinite(photons(:))));
end

function testVernierOisCreate(testCase)
vparams(2) = vernierP;
vparams(2).name = 'offset';
vparams(2).bgColor = 0;
vparams(1) = vparams(2);
vparams(1).barWidth = 0;
vparams(1).bgColor = 0.5;
sparams.fov = 0.5;
stimWeights = ones(1, 5) / 5;

ois = oisCreate('vernier', 'add', stimWeights, ...
    'testParameters', vparams, 'sceneParameters', sparams);

testCase.verifyClass(ois, 'oiSequence');
testCase.verifyEqual(ois.length, numel(stimWeights));
end

function testInvalidTypeErrors(testCase)
% An oisType not in {'harmonic','vernier','impulse'} must fail at parse time.
stimWeights = ones(1, 5);
threw = false;
try
    oisCreate('notavalidtype', 'blend', stimWeights);
catch
    threw = true;
end
testCase.verifyTrue(threw, 'Expected error for invalid oisType.');
end
