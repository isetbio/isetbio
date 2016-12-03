%% t_rgcImpulseResponsePillow
%
% Plot the pillow temporal impulse response function used in the
% conventional coupled-GLM.
%
% As implemented, the Pillow impulse response function produces a fixed
% shape, and the time axis shifts depending on the filter duration.  It
% could be written to produce the same curve.  Not sure what is intended.
%
% BW ISETBIO Team, 2016

%%

clear params

params.filterDuration = 0.2; params.samplingTime = 0.002;
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

vcNewGraphWin; 
plot(timeAxis,rgcFilter,'r-x'); xlabel('Sec'); grid on

params.filterDuration = 0.3; params.samplingTime = 0.005;
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

hold on; 
plot(timeAxis,rgcFilter,'b-o'); xlabel('Sec'); grid on

%%