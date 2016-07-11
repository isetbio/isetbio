%% Cone mosaic object - introduction
%
%

%%
ieInit

%% Build a scene and oi for computing

% s = sceneCreate('vernier');
% s.distance = 1;

s = sceneCreate('rings rays');
vcAddObject(s);
% s = sceneCreate('slanted bar');
% fname = fullfile(isetRootPath,'data','images','rgb','eagle.jpg');
% s = sceneFromFile(fname,'rgb');

s = sceneSet(s,'fov',2);

oi = oiCreate;
oi = oiCompute(oi,s);
vcAddObject(oi); % oiWindow;

%% Build a default cone mosaic and compute the OI

cMosaic = coneMosaic;  % Create the object
% cMosaic.rows = 100; cMosaic.cols = 120;
cMosaic.rows = 144; cMosaic.cols = 176;
cMosaic.emGenSequence(500);
cMosaic.compute(oi,'currentFlag',true);

% @BW/@HJ
% There appears to be a mismatch between the oi field of view and the size
% of the cone mosaic.  The cone mosaic is smaller than the field of view
% set in the scene and oi.


% Show the window
% cMosaic.guiWindow;

% Examine the outer segment current
% cMosaic.plot('movie absorptions','vname','deleteme.avi','step',5);

%% To compute the bipolar response
bp = bipolar(cMosaic.os);
bp.compute(cMosaic.os);

% bp.plot('response');
% bp.plot('movie response');

%% To compute an RGC response
% Build rgc
clear params
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 1;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

%
ir = irCreate(bp, params);
ir.mosaicCreate('model','GLM','type','on midget');
% Some testing of the new rgc get and set organization
% ir.mosaic{1}.get('cell type')
% ir.mosaic{1}.get('rf diameter')
% ir.mosaic{1}.get('srf center')
% ir.mosaic{1}.get('srf surround')
% ir.mosaic{1}.get('generator function')
%
% ir.mosaic{1}.set('number trials',1)

% Number of repeated trials
% ir.mosaicCreate('model','lnp','type','on midget');
% ir.mosaic{1}.set('numberTrials',3);

fprintf('Cell array size: %d x %d\n',ir.mosaic{1}.get('mosaic size'));

%% Compute RGC response
tic;
ir = irCompute(ir, bp, 'coupling',false);
toc

lastTime = ir.mosaic{1}.get('last spike time');

ir.mosaic{1}.set('dt',1);
psth = ir.mosaic{1}.get('psth');

clear params
params.vname = fullfile(isetbioRootPath,'local','vernier.avi'); 
param.FrameRate = 5; params.step = 2; params.show = false;
ieMovie(psth,params);

%%
% irPlot(ir, 'mosaic');
% irPlot(ir, 'linear');
irPlot(ir, 'raster','cell',[5,7]);
irPlot(ir, 'raster');   % Can be very long and painful.  Fix it in the plot
irPlot(ir, 'mosaic'); 

% irPlot(ir, 'psth');

%%