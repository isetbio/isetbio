function [settings, primaries, predictedSpds] = SpdToSettings(calOrCalStruct, targetSpds, varargin)
% Converts spectra into settings.
%
% Syntax:
%    [settings, primaries] = SpdToSettings(calOrCalStruct, targetSpds)
%    [settings, primaries] = SpdToSettings(calOrCalStruct, targetSpds, 'lambda', '001')
%
% Description:
%    Essentially a convience wrapper around SpdToPrimary and
%    PrimaryToSettings.  Returns settings, primaries and predicted
%    spectra.
%
% Input:
%    calOrCalStruct               - Calibration struct or object.
%    targetSpds (nWls x nSpectra) - Target spectra.  Should be on the same wavelength
%                                   spacing and power units as the PR-650 field of the
%                                   calibration structure.
%
% Output:
%    settings (nPrimaries x nSpectra) - The [0,1], gamma corrected power level for each
%                                       effective primary of the OneLight.  Each column is
%                                       a single set of primary values.
%    primaries (nPrimaries x nSpectra) - The primary values.
%    predictedSpds (nWls x nSpectra) -  The predicted spectra 
%
% Optional key/value pairs:
%   'lambda' - scalar (default 0.1) - Determines how much smoothing we apply to the settings.
%   'verbose' - true/false (default false) - Enables/disables verbose diagnostic information.

% 09/11/21  dhb  Wrote from OLSpdToSettings

%% Parse the input
p = inputParser;
p.addParameter('verbose', false, @islogical);
p.addParameter('lambda', 0.1, @isscalar);
p.parse(varargin{:});
params = p.Results;

%% Make sure we have @CalStruct object that will handle all access to the calibration data.
%
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
    error('The input (calOrCalStruct) is not a cal struct.');
end

%% Check wavelength sampling
ambient = scalStructOBJ.get('P_ambient');

S = calStructOBJ.get('S');
wls = SToWls(S);
nWls = size(targetSpds,1);
if (nWls ~= S(3))
    error('Wavelength sampling inconsistency between passed spectrum and calibration');
end

% Convert the spectra into primaries.
nPrimaries = size(calOrCalStruct.computed.pr650M,2);
numSpectra = size(targetSpds, 2);
primaries = zeros(nPrimaries, numSpectra);
for i = 1:numSpectra
	% Convert to primaries.
	primaries(:,i) = SpdToPrimary(calOrCalStruct, targetSpds(:,i), 'lambda',params.lambda, 'verbose', params.verbose);
    
    % Convert to spectra
    predictedSpds(:,i) = PrimaryToSpd(calOrCalStruct,primaries(:,i));
end

% Convert from primaries to settings.
settings = PrimaryToSettings(calOrCalStruct, primaries);
