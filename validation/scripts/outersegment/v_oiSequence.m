function varargout = v_oiSequence(varargin)
%
% Validate simulations using the oiSequence method
%  
%  * Create a small oiSequence
%  * Check the coneMosaic
%  * Adjust the time base of the oiSequence and repeat
%
% The purpose is to show that we can vary the oiSequence time base and the
% cMosaic integration time and have things roughly work out.
%
% This script test oiSequence and computeForOISequence() as well as their
% dependent functions.
%
% NPC, ISETBIO TEAM, 2016

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

% Reproduce identical random number
rng(1, 'combRecursive');

%% Produces the harmonic.  Need to lock this down as a function.

stimWeights = ieScale(fspecial('gaussian',[1,10],3),0,1);
weights = [zeros(1, 5), stimWeights, zeros(1, 5)];
 
hparams(1) = harmonicP;
hparams(2) = hparams(1); hparams(2).contrast = 0;
sparams.fov = 0.3;
ois = oisCreate('harmonic','blend',weights, 'testParameters',hparams,'sceneParameters',sparams);
% ois.visualize('movie illuminance');

% Assert values from December 22, 2017 (DHB)
% Updated by hand from previous values because of changes to optics code.
tolerance = 1E-3;

% Check the oiCreate part
photonsFixed = oiGet(ois.oiFixed,'photons');
quantityOfInterest = sum(photonsFixed(:))/1.3085e+19 - 1;
UnitTest.assertIsZero(quantityOfInterest,'oiFixed photons',tolerance);

photonsModulated = oiGet(ois.oiModulated,'photons');
quantityOfInterest = sum(photonsModulated(:))/1.3085e+19 - 1;
UnitTest.assertIsZero(quantityOfInterest,'oiModulated photons',tolerance);

% This tests the generation of the sequence because the sequence is built
% by coneMosaic.compute
cMosaic = coneMosaic;
cMosaic.noiseFlag = 'frozen';
cMosaic.os.noiseFlag = 'frozen';
cMosaic.setSizeToFOV(0.2);
cMosaic.integrationTime = ois.timeAxis(2);  % This is the integration time
cMosaic.emGenSequence(length(ois.timeAxis),'rseed',1);
cMosaic.compute(ois);

tolerance = 2E-2;
totalAbsorptions = sum(cMosaic.absorptions(:));

% Changed from the value 116255 when updating the eye movement model
% to fixationalEM in emGenSequence 
quantityOfInterest = (totalAbsorptions/127714) - 1;   
UnitTest.assertIsZero(quantityOfInterest,'coneMosaic absorptions',tolerance);

% assert((sum(cMosaic.absorptions(:))/121450) - 1 < 1e-3);

% Unit test validation data
UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_oiSequence *****');
UnitTest.validationData('photonsFixed', photonsFixed);
UnitTest.validationData('photonsModulated', photonsModulated);
UnitTest.validationData('totalAbsorptions', totalAbsorptions, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'totalAbsorptions', 600);

end


