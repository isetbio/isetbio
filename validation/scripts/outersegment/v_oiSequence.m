function varargout = v_oiSequence(varargin)
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
rng('default'); rng(1);

%% Produces the harmonic.  Need to lock this down as a function.

stimWeights = ieScale(fspecial('gaussian',[1,10],3),0,1);
weights = [zeros(1, 5), stimWeights, zeros(1, 5)];
 
hparams(1) = harmonicP;
hparams(2) = hparams(1); hparams(2).contrast = 0;
sparams.fov = 0.3;
ois = oisCreate('harmonic','blend',weights, 'tparams',hparams,'sparams',sparams);
% ois.visualize;

% Assert values from November 19, 2016 (BW)
tolerance = 1E-3;

% Check the oiCreate part
photons = oiGet(ois.oiFixed,'photons');
quantityOfInterest = sum(photons(:))/1.3379e+19 - 1;
UnitTest.assertIsZero(quantityOfInterest,'oiFixed photons',tolerance);

photons = oiGet(ois.oiModulated,'photons');
quantityOfInterest = sum(photons(:))/1.3379e+19 - 1;
UnitTest.assertIsZero(quantityOfInterest,'oiModulated photons',tolerance);

% This tests the generation of the sequence because the sequence is built
% by coneMosaic.compute
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.2);
cMosaic.integrationTime = ois.timeAxis(2);  % This is the integration time
cMosaic.emGenSequence(length(ois.timeAxis));
cMosaic.compute(ois);

tolerance = 1E-2;
quantityOfInterest = (sum(cMosaic.absorptions(:))/121450) - 1;
UnitTest.assertIsZero(quantityOfInterest,'coneMosaic abosprtions',tolerance);

% assert((sum(cMosaic.absorptions(:))/121450) - 1 < 1e-3);

end


