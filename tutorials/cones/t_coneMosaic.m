%% Cone mosaic object - introduction
%
%

%%
ieInit

%% Build a scene and oi for computing

% s = sceneCreate('rings rays');
s = sceneCreate('slanted bar');
% fname = fullfile(isetRootPath,'data','images','rgb','eagle.jpg');
% s = sceneFromFile(fname,'rgb');

s = sceneSet(s,'fov',2);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); oiWindow;

%% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic;                     % Create the object
cMosaic.emGenSequence(500);
cMosaic.compute(oi,'currentFlag',true);   % The current is computed by default anyway

% Show the window
cMosaic.guiWindow;                      % Note: Default image should be mean absorptions

% Examine the outer segment current
% cMosaic.plot('current time series');

%% To compute the bipolar response

bp = bipolar(cMosaic.os);

bp.compute(cMosaic.os);

bp.plot('response');
%% To compute an RGC response

% Build rgc

clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

%
ir = irCreate(bp, params);
ir.mosaicCreate('model','glm','type','off parasol');

% Compute RGC response
ir = irCompute(ir, bp);

%%
% irPlot(ir, 'mosaic');
% irPlot(ir, 'linear');
irPlot(ir, 'raster');
% irPlot(ir, 'psth');

%%