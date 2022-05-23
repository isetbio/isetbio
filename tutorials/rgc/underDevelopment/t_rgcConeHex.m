% t_rgcConeMosaicHex
% 
% Compute an RGC mosaic response from the coneMosaicHex object for a moving
% bar stimulus.
% 
% This tutorial generates a moving bar stimulus and computes a
% coneMosaicHex response. The cone mosaic is passed on to a bipolar mosaic
% and then an RGC mosaic. The computed responses are visualized along the
% way. This a demonstration that the bipolar and RGC objects can accept and
% compute with the coneMosaicHex object.
% 
% 9/2016 JRG/BW (c) isetbio team

%% Build moving bar stimulus, oi and coneMosaicHex
clear
showMovieFlag = 1;
params.barWidth = 10;
params.fov      = 0.3;
% params.os = 'biophys';
params.os = 'hex';

% Alternatively, we will be able to do this:
%
%   rd = RdtClient('isetbio');
%   rd.crp('/resources/data/istim');
%   rd.listArtifacts('print',true);
%  And then upload the relevant one
%

iStim = ieStimulusBar(params);  % Full params are returned in iStim

%% Visualize coneMosaicHex responses

% iStim.cMosaic.visualizeGrid();

isomerizationsBar = iStim.cm.absorptions;
% iStim.cMosaic.visualizeActivationMaps(...
%     isomerizationsBar(:,:,1), ...                                  % the response matrix
%        'mapType', 'modulated hexagons', ...                          % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
%     'signalName', 'isomerizations (R*/cone/integration time)', ...   % colormap title (signal name and units)
%       'colorMap', jet(1024), ...                                     % colormap to use for displaying activation level
%     'figureSize', [1550 950] ...                                     % figure size in pixels
%     );

% iStim.cMosaic.window;

%% Compute bipolar response

bpParams.cellType = 'offmidget';
bp = bipolar(iStim.cm, bpParams);
bp.set('sRFcenter',10);
bp.set('sRFsurround',0);
bp.compute(iStim.cm);
bp.plot('movie response')

%% Compute RGC response
clear params innerRetinaSU
cellType = 'offMidget';
% cellType = 'offParasol';
params.name = 'macaque phys';
params.eyeSide = 'left';
params.eyeRadius = 0;
params.eyeAngle = 0; ntrials = 0;

% Create RGC object
innerRetinaSU = ir(bp, params);
innerRetinaSU.mosaicCreate('type',cellType,'model','GLM');

nTrials = 1; innerRetinaSU = irSet(innerRetinaSU,'numberTrials',nTrials);

% Compute the inner retina response
innerRetinaSU = irCompute(innerRetinaSU, bp);

%% Show output

% Make the PSTH movie
lastTime = innerRetinaSU.mosaic{1}.get('last spike time');
innerRetinaSU.mosaic{1}.set('dt',1);
psthTest = innerRetinaSU.mosaic{1}.get('psth');


% View movie of RGC linear response
clear vParams
vParams.FrameRate = 5; vParams.step = 2; vParams.show = true;

vcNewGraphWin;
lResponse = innerRetinaSU.mosaic{1}.get('response linear');
% ieMovie(lResponse,vParams);

% View movie of PSTH for mosaic
if showMovieFlag
    steadyStateFrame = 50; % Get rid of transient spiking
    vcNewGraphWin;
    ieMovie(psthTest(:,:,steadyStateFrame:end),vParams);
end

%%
