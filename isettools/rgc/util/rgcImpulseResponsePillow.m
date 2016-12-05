function [rgcFilter,timeAxis]  = rgcImpulseResponsePillow(varargin)
% Build the temporal impulse response used by Pillow et l.
% 
%    [rgcFilter, timeAxis] = ieRGCTIRPillow([params])
%
% At present, the Pillow impulse response function produces a fixed shape,
% but the time axis scales depending on the filter duration.  It could be
% written to produce the same curve.  Not sure what is intended.
%
% The default parameters in the code from Pillow has filterLength (which
% means duration) of 200ms.
%
% Inputs
%  filterDuration - Duration in seconds
%  samplingTime - Duration of each sample bin in seconds
%
% Returns:
%   rgcFilter - The Pillow filter
%   timeAxis  - Sample times in seconds
%
% Example
%    params.filterDuration = 0.2; params.samplingTime = 0.002;
%    [rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);
%    vcNewGraphWin; plot(timeAxis,rgcFilter); xlabel('Sec'); grid on
%
%    params.filterDuration = 0.3; params.samplingTime = 0.005; 
%    [rgcFilter,timeAxis]  = rgcImpulseResponsePillow(params);
%    vcNewGraphWin; plot(timeAxis,rgcFilter,'-o'); xlabel('Sec'); grid on
%
% BW, ISETBIO Team, 2016

%%
p = inputParser;

p.addParameter('filterDuration',0.2,@isnumeric);    % Duration in sec
p.addParameter('samplingTime',0.002,@isnumeric);  % Sample times in sec

p.parse(varargin{:});
filterDuration = p.Results.filterDuration;
samplingTime = p.Results.samplingTime;

%% Compute the curve, respecting temporal sample

nkt = ceil(filterDuration/samplingTime);  % Number of time bins

% tk = timeAxis;
tk = (0:nkt-1)';
timeAxis = tk*samplingTime;

b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf

rgcFilter = (k1 - k2./1.5);
rgcFilter = 1.2*(rgcFilter./max(rgcFilter));

end

%% Original Pillow code from buildTemporalImpulseResponse.m
%
% samplingTime = 1;
% DTsim = .01; % Bin size for simulating model & computing likelihood.
% nkt = 20;  % Number of time bins in filter;
% DTsim = samplingTime; % Bin size for simulating model & computing likelihood.
% filterLength = 0.2;
% nkt = ceil(filterLength/samplingTime);  % Number of time bins in filter;
% timeAxis = samplingTime:samplingTime:filterLength;
% tk = [0:nkt-1]';
% b1 = nkt/32; b2 = nkt/16;
% k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
% k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
% k = (k1-k2./1.5);
% k = 1.2*(k./max(k));
% plot(k)
%