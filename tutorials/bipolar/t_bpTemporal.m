%% t_bpTemporal
%
% Illustrating the reasons for choosing the bipolar filters.
%
% In bipolarFilter we get the Pillow impulse response function.
%
% BW ISETBIO Team, 2016

%%
ieInit
% clear variables

%%  First illustrate a small cone mosaic
rows = 5; cols = 5; tSamples = 100;
cMosaic = coneMosaic('size',[rows,cols]);

cMosaic.absorptions = single(randi(10,rows,cols,tSamples));
cMosaic.integrationTime = 0.002;
cMosaic.emGenSequence(tSamples);
cMosaic.window;

%%  Compute the os linear filters and then make the bipolar

cMosaic.os.linearFilters(cMosaic);

bp = bipolarCreate(cMosaic);

% The bipolarFilter routine tries to create a filter so that os convolved
% with bipolar matches the Pillow filter.
bpFilter = bipolarFilter(bp, cMosaic,'graph',true);

%%
% If you used a cMosaic with a real time varying signal, then you can do
% this to see the bipolar filter at the cone mosaic sample.
vcNewGraphWin; plot(cMosaic.timeAxis,bpFilter,'o-');

%%