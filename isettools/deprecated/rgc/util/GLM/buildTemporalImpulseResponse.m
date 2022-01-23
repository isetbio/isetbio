function [k, timeAxis] = buildTemporalImpulseResponse(samplingTime)
% Builds the temporal filter function for RGCs in an rgcMosaic.
%
% Syntax:
%   ir.tCenter{rgbInd, 1} = buildTemporalImpulseResponse(samplingTime);
%
% Description:
%    Create a default (temporal) stimulus filter using the gamma PDF.
%
%    This function contains examples of usage inline. To access these, type
%    'edit buildTemporalImpulseResponse.m' into the Command Window.
%
% Inputs:
%    samplingTime - Numeric. The sampling time from the outersegment object
%
% Outputs:
%    k            - Array. A Nx1 array where N = 0.2/samplingTime. This
%                   array contains the temporal impulse response.
%    timeAxis     - Array. A 1xN array containing the relevant time axis
%                   for the temporal impulse response.
%
% Optional key/value pairs:
%    None.
%
% References:
%    * This function incorporates code by Pillow available under the GNU
%      General Public License, at:
%      http://pillowlab.princeton.edu/code_GLM.html
%

% History:
%    09/XX/15  JRG  Created
%    06/10/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define obj before it could
    % possibly work.
    obj.tCenter{rgbInd, 1} = multFactor * ...
        buildTemporalImpulseResponse(samplingTime);
%}

%% Build temporal impulse response
% From code by J. Pillow:
% 
% samplingTime = 1;
% DTsim = .01; % Bin size for simulating model & computing likelihood.
% nkt = 20;  % Number of time bins in filter;
% % Bin size for simulating model & computing likelihood (seconds)
% DTsim = samplingTime;

filterLength = 0.2;   % Seconds
nkt = ceil(filterLength / samplingTime);  % Number of time bins in filter;
timeAxis = samplingTime:samplingTime:filterLength;

tk = (0:nkt - 1)';
b1 = nkt / 32;
b2 = nkt / 16;
k1 = 1 / (gamma(6) * b1) * (tk / b1) .^ 5 .* exp(-tk / b1);  % Gamma pdfn
k2 = 1 / (gamma(6) * b2) * (tk / b2) .^ 5 .* exp(-tk / b2);  % Gamma pdf

k = (k1 - k2 ./ 1.5);
k = 1.2 * (k ./ max(k));

% plot(k)
end