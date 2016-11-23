function varargout = v_cmosaic(varargin)
%
% Simple cone mosaic calculation.  We will systematically change parameters
% and see that the results are stable. 
%
% BW, ISETBIO Team Copyright 2016

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Create appropriate structures
scene = sceneCreate('uniform ee');
scene = sceneSet(scene,'fov',1);
scene = sceneAdjustLuminance(scene,200);
oi = oiCreate('human');
oi = oiCompute(oi,scene);

%%
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(0.5*sceneGet(scene,'fov'));
cMosaic.compute(oi);
% cMosaic.window;

sumA = sum(cMosaic.absorptions(:));

% How I calculated
% A = 0;
% for ii=1:10; cMosaic.compute(oi); A = A + sum(cMosaic.absorptions(:)); end
% v = A/10;

v = 3.7358e+06;
UnitTest.assert((sumA < (v + 4*sqrt(v))),'Not too hot... ');
UnitTest.assert((sumA > (v - 4*sqrt(v))),'Not too cold.. ');

%% Plot
if (runTimeParams.generatePlots)
    
end

%% Do the whole thing again, but with the macular pigment set to zero
cMosaic.macular.density=0;
cMosaic.compute(oi);
% cMosaic.window;

sumA = sum(cMosaic.absorptions(:));
% A = 0;
% for ii=1:10; cMosaic.compute(oi); A = A + sum(cMosaic.absorptions(:)); end
% v = A/10;

v = 4.1820e+06;
UnitTest.assert((sumA < (v + 4*sqrt(v))),'Not too many.. ');
UnitTest.assert((sumA > (v - 4*sqrt(v))),'Not too few... ');
%% Plot
if (runTimeParams.generatePlots)
    
end

end

