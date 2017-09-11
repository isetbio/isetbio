%% Spike triggered average estimates of RGC
%
% Not tested or checked
% Needs more comments.
%
% Check that the grid appears as a grid on the RGC response.
%
% JRG/BW, ISETBIO Team, 2017

%%
ieInit;

%% Scene, oi, grid, cone mosaic

imSize = 32; 

% Should be long when not debugging.  
% Otherwise STA doesn't look good.
tSteps = 500;    
coneParams.fov = 0.1;
wnStimulus = 255*round(rand(imSize,imSize,tSteps));

% Stimulus looks off.  Saturated or something.
tic
ieStim = ieStimulusMovieCMosaic(wnStimulus,coneParams);
cMosaic = ieStim.cMosaic;
toc
%% Make the biplar layer with just one mosaic 

clear bpL bpMosaicParams
bpL = bipolarLayer(cMosaic);

ii = 1;
bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{ii} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);
bpL.mosaic{ii}.compute;

% bpL.window;
%% Make the RGC layer and show it

clear rgcL rgcParams
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 2;

% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget');
rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget',rgcParams,'bipolarScale',100,'bipolarContrast',1);
rgcL.compute;

% rgcL.window;

%%  Set the mean of every stimulus to zero

% If there is a spike 

spikeResp = mosaicSpikes(rgcL);
nTimes = 15;
% If we make make this a vector, it will do it for all multiple cells.
cellList = 10;

% The stimulus is in XW format (really XT).  So whenever there is a spike
% we add that stimulus in.
[wnXW,imRows,imCols] = RGB2XWFormat(wnStimulus);

% Cell identity, image size, nTime steps
sta = zeros(length(cellList),size(wnXW,1),nTimes);

% The zero mean stimulus
stimzm = (single(wnXW)-(ones(size(wnXW,1),1)*mean(wnXW,1)));

% This could be simpler.  It is a matrix multiplication, really.  The
% timing is a little odd.  Re-write to make the math clearer.
for j = 1:nTimes   % Steps back in time
    sta(:,:,j) = ((stimzm(:,1:end-(j-1))))*(spikeResp(cellList,j:size(stimzm,2))'); 
end

staMovie = XW2RGBFormat(squeeze(sta(1,:,:)),imRows,imCols);

%% Not correct at this point.  Figure out why.

ieMovie(staMovie);

%%
