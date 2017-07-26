%% t_rgcImpulseResponsePillow
%
% Plot the pillow temporal impulse response function used in the
% conventional coupled-GLM.
%
% As implemented, the Pillow impulse response function produces a fixed
% shape, but the time axis scales depending on the filter duration.  
%
% It could be written to produce the same curve.  Not sure what is
% intended. 
%
% BW ISETBIO Team, 2016

%%

clear params

params.filterDuration = 0.2; params.samplingTime = 0.001;
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

vcNewGraphWin; 
plot(timeAxis,rgcFilter,'r-x'); xlabel('Sec'); grid on

params.filterDuration = 0.2; params.samplingTime = 0.005;
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

hold on; 
plot(timeAxis,rgcFilter,'b-o'); xlabel('Sec'); grid on
legend('Sampling 1 ms','Sampling 5 ms');

%%

clear params

params.filterDuration = 0.1; params.samplingTime = 0.003;
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

vcNewGraphWin; 
plot(timeAxis,rgcFilter,'r-x'); xlabel('Sec'); grid on

params.filterDuration = 0.2; 
[rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);

hold on; 
plot(timeAxis,rgcFilter,'b-o'); xlabel('Sec'); grid on
legend('Duration 100ms','Duration 200ms');

%%