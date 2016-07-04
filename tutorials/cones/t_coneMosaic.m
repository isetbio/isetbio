%% Cone mosaic object - introduction
%
%

%%
ieInit

%% Build a scene and oi for computing

s = sceneCreate('vernier');
s.distance = 1;

% s = sceneCreate('rings rays');
% s = sceneCreate('slanted bar');
% fname = fullfile(isetRootPath,'data','images','rgb','eagle.jpg');
% s = sceneFromFile(fname,'rgb');

s = sceneSet(s,'fov',2);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); oiWindow;

%% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic;                     % Create the object
% cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(500);
cMosaic.compute(oi,'currentFlag',true);   % The current is computed by default anyway

% Show the window
cMosaic.guiWindow;                      % Note: Default image should be mean absorptions

% Examine the outer segment current
% cMosaic.plot('current time series');

%% To compute the bipolar response

bp = bipolar(cMosaic.os);

bp.compute(cMosaic.os);

% bp.plot('response');
clear params
params.vname = tempname; param.FrameRate = 5; params.step = 2; params.show = true;
bp.plot('movie response',params);

%% To compute an RGC response

% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 0.5;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

%
ir = irCreate(bp, params);
ir.mosaicCreate('model','lnp','type','on midget');

% Number of repeated trials
ir.mosaic{1}.set('numberTrials',3);

fprintf('Cell array size: %d x %d\n',ir.mosaic{1}.get('mosaicsize'));
% Compute RGC response
ir = irCompute(ir, bp);
lastTime = ir.mosaic{1}.get('last spike time');

psth = ir.mosaic{1}.get('psth','dt',1);

clear params
params.vname = 'vernier'; param.FrameRate = 5; params.step = 2; params.show = true;
ieMovie(psth,params);

%%
% irPlot(ir, 'mosaic');
% irPlot(ir, 'linear');
irPlot(ir, 'raster','cell',[5,7]);
irPlot(ir, 'raster');   % Can be very long and painful.  Fix it in the plot

% irPlot(ir, 'psth');

%%