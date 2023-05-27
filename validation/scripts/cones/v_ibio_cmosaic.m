function varargout = v_cmosaic(varargin)
%
% Simple rectangular cone mosaic calculation.
%
% We will systematically change parameters and see that the results are stable. 
%
% BW, ISETBIO Team Copyright 2016

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)


%% Initialize
ieInit;

%% Reproduce identical random numbers
rng('default'); rng(1);

%% Create appropriate structures
scene = sceneCreate('uniform ee');
scene = sceneSet(scene,'fov',1);
scene = sceneAdjustLuminance(scene,200);
oi = oiCreate('human');
oi = oiCompute(oi,scene);

%% Creat the mosaic and compute isomerizations (which are sometimes called absorptions)
% This doesn't use the default integration time.

cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.5*sceneGet(scene,'fov'));
cMosaic.integrationTime = 0.050;    
cMosaic.compute(oi);

%% Check that they are within Poisson noise of what we expected at some point when this was created.
sumA = sum(cMosaic.absorptions(:));
expectedSumA = 3738338;
UnitTest.assert((sumA < (expectedSumA + 4*sqrt(expectedSumA))),'Not too hot... ');
UnitTest.assert((sumA > (expectedSumA - 4*sqrt(expectedSumA))),'Not too cold.. ');

%% Do the whole thing again, but with the macular pigment set to zero
cMosaic.macular.density = 0;
cMosaic.compute(oi);
sumA = sum(cMosaic.absorptions(:));
expectedSumA = 4.1820e+06;
UnitTest.assert((sumA < (expectedSumA + 4*sqrt(expectedSumA))),'Not too many.. ');
UnitTest.assert((sumA > (expectedSumA - 4*sqrt(expectedSumA))),'Not too few... ');

%% Plot if there are any
if (runTimeParams.generatePlots)
    
end

end

