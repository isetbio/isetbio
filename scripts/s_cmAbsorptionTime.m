%% Check absorption shapes for different cone mosaic sizes
%
%    s_cmAbsorptionSize
%
% Check that the dimensions of the absorptions are correct for different
% size cone arrays and temporal samples.  We test 2D and 1D spatial arrays
% with a single time point and with multiple time points.
%
% BW, ISETBIO Team, 2016

%% 
ieInit

scene = sceneCreate('uniform ee');
oi = oiCreate;
oi = oiCompute(oi,scene);

cMosaic = coneMosaic('size',[10 10]);
cMosaic.integrationTime = 0.05;
cMosaic.emGenSequence(30);
cMosaic.compute(oi);

% In this case, we have a matrix.  size(
cMosaic.window;

%%
cMosaic.integrationTime = 0.01;
cMosaic.emGenSequence(30);
cMosaic.compute(oi);
cMosaic.window;



%% One dimensional cone array, 10 times

cMosaic = coneMosaic('size',[1 100]);
cMosaic.emGenSequence(10);
cMosaic.compute(oi);
cMosaic.window;

assert(isequal(size(cMosaic.absorptions),[1,100,10]));

%% Two dimensional cone array, 20 times

cMosaic = coneMosaic('size',[20 20]);
cMosaic.emGenSequence(20);
cMosaic.compute(oi);
cMosaic.window;

assert(isequal(size(cMosaic.absorptions),[20,20,20]));

%%

