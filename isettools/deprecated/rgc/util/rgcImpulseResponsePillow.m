function [rgcFilter, timeAxis] = rgcImpulseResponsePillow(varargin)
% Build the temporal impulse response used by Pillow et l.
%
% Syntax:
%   [rgcFilter, timeAxis] = rgcImpulseResponsePillow([varargin])
%
% Description:
%    At present, the Pillow impulse response function produces a fixed
%    shape, but the time axis scales depending on the filter duration. It
%    could be written to produce the same curve.  Not sure what is
%    intended.
%
%    The default parameters in the code from Pillow has filterLength (which
%    means duration) of 200ms.
%
%    This function contains examples of usage inline. To access these, type
%    'edit rgcImpulseResponsePillow.m' into the Command Window.
%
% Inputs:
%    None required.
%
% Outputs:
%    rgcFilter      - The Pillow filter
%    timeAxis       - Sample times in seconds
%
% Optional key/value pairs:
%    filterDuration - Numeric. The duration in seconds. Default 0.2.
%    samplingTime   - Numeric. The duration of each sample bin in seconds.
%                     Default 0.002.
%    cellType       - String. The parasol cell type. Default 'onParasol'.
%                     Options include the following: onparasol, offparasol,
%                     ondiffuse, offdiffuse, sbc, onsbc, and then
%                     onmidget, and offmidget.
%

% History:
%    XX/XX/16  BW   ISETBIO Team, 2016
%    06/05/19  JNM  Documentation pass

% Example
%{
    params.filterDuration = 0.2;
    params.samplingTime = 0.002;
    [rgcFilter, timeAxis] = rgcImpulseResponsePillow(params);
    vcNewGraphWin;
    plot(timeAxis, rgcFilter);
    xlabel('Sec');
    grid on;

    params.filterDuration = 0.3;
    params.samplingTime = 0.005; 
    [rgcFilter, timeAxis] = rgcImpulseResponsePillow(params);
    vcNewGraphWin;
    plot(timeAxis, rgcFilter, '-o'); 
    xlabel('Sec');
    grid on;
%}

%%
p = inputParser;
p.addParameter('filterDuration', 0.2, @isnumeric);  % Duration in sec
p.addParameter('samplingTime', 0.002, @isnumeric);  % Sample times in sec
p.addParameter('cellType', 'onparasol', @ischar);   % Sample times in sec
p.parse(varargin{:});
filterDuration = p.Results.filterDuration;
samplingTime = p.Results.samplingTime;
cellType = p.Results.cellType;
%% Compute the curve, respecting temporal sample
nkt = ceil(filterDuration/samplingTime);  % Number of time bins

% tk = timeAxis;
tk = (0:nkt - 1)';
timeAxis = tk*samplingTime;

switch cellType
    case {'onparasol', 'offparasol', 'ondiffuse', 'offdiffuse', ...
            'sbc', 'onsbc'}
        b1 = nkt / 32;
        b2 = nkt / 16;
        c1 = 1;
        c2 = 1 / 1.5;        

        % Gamma pdfn and gamma pdf (1 and 2 below respectively)
        k1 = 1 / (gamma(6) * b1) * (tk / b1) .^ 5 .* exp(-tk / b1);
        k2 = 1 / (gamma(6) * b2) * (tk / b2) .^ 5 .* exp(-tk / b2);
        rgcFilter = (c1 * k1 - c2 * k2);
        rgcFilter = 1.2 * (rgcFilter ./ max(rgcFilter));
    case {'onmidget', 'offmidget'}
        % fit from apricot data set
        tk = tk / 100; % need to shorten time base for good gamma fit
        % b1 = 0.08471;
        % b2 = 0.3827;  % at nkt = 400
        b1 = nkt / (400 / 0.08471);
        b2 = nkt / (400 / 0.3827);
        c1 = 1.199;
        c2 = 0.3188;        

        % Gamma pdfn and gamma pdf (1 and 2 below respectively)
        k1 = 1 / (gamma(6) * b1) * (tk / b1) .^ 5 .* exp(-tk / b1);
        k2 = 1 / (gamma(6) * b2) * (tk / b2) .^ 5 .* exp(-tk / b2);
        rgcFilter = (c1 * k1 - c2 * k2);
        rgcFilter = 1.2 * (rgcFilter ./ max(rgcFilter));
end
end

%% Original Pillow code from buildTemporalImpulseResponse.m
%{
    samplingTime = 1;
    DTsim = .01; % Bin size for simulating model & computing likelihood.
    nkt = 20;    % Number of time bins in filter;
    % Bin size for simulating model & computing likelihood.
    DTsim = samplingTime;
    filterLength = 0.2;
    nkt = ceil(filterLength/samplingTime);  % Number of time bins in filter
    timeAxis = samplingTime:samplingTime:filterLength;
    tk = [0:nkt - 1]';
    b1 = nkt / 32;
    b2 = nkt / 16;
    % Gamma pdfn and gamma pdf (1 and 2 below respectively)
    k1 = 1 / (gamma(6) * b1) * (tk / b1) .^ 5 .* exp(-tk / b1);
    k2 = 1 / (gamma(6) * b2) * (tk / b2) .^ 5 .* exp(-tk / b2);
    k = (k1 - k2 ./ 1.5);
    k = 1.2 * (k ./ max(k));
    plot(k)
%}