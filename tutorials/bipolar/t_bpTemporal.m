%% t_bpTemporal
%
%

rows = 5; cols = 5; tSamples = 100;
cMosaic = coneMosaic('size',[rows,cols]);

cMosaic.absorptions = single(randi(10,rows,cols,tSamples));
cMosaic.integrationTime = 0.002;
cMosaic.emGenSequence(tSamples);
cMosaic.window;

cMosaic.os.linearFilters(cMosaic);
bp = bipolarCreate(cMosaic);

% The bipolarFilter routine tries to create a filter so that os convolved
% with bipolar matches the Pillow filter.
bpFilter = bipolarFilter(bp, cMosaic,'graph',true);

% If you used a cMosaic with a real time varying signal, then you can do
% this.  See t_bpTemporal
vcNewGraphWin; plot(cMosaic.timeAxis,bpFilter,'o-');

%%