function idx = ieFindWaveIndex(wave, waveVal, perfect)
% Returns a (0/1) vector of indices s.t. wave(idx) matches a waveVal entry.
%
% Syntax:
%   idx  = ieFindWaveIndex(wave, waveVal, [perfect])
%
% Description:
%    We want to address only some wavebands in a radiance data set (e.g., 
%    photons) we find the relevant indices in wave by this call. 
%
%    If perfect = 1, this routine uses the Matlab function ismember(). 
%
%    If perfect = 0, we accept a closest match, say we want the closest
%    value. In this case run, the same wave valuel may match two waveVal
%    entries, and there will be different vector lengths returned. We
%    announce this mis-match condition.
%
% Inputs:
%    wave    - The wavelengths
%    waveVal - Desired 
%    perfect - (Optional) Boolean indicating whether or not to look for
%              perfect value match(es) for waveVal in wave. Default is 1. 
%
% Outputs:
%    idx     - A vector of booleans showing which of the instances of wave
%              match a waveVal entry
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/30/17  jnm  Formatting & fix example
%

% Examples:
%{
    scene = sceneCreate;
    wave = sceneGet(scene, 'wave');
    waveVal = [500, 600];
    idx = ieFindWaveIndex(wave, waveVal);
%}

if notDefined('wave'), error('Must define list of all wavelengths'); end
if notDefined('waveVal'), error('Must define wavelength values'); end
if notDefined('perfect'), perfect = 1; end

if perfect
    % Find only perfect matches
    idx = logical(ismember(wave, waveVal));
else
    idx = false(1, length(wave));   % Assume not a member
    % For each waveVal, find the index in wave that is closest.
    for ii=1:length(waveVal)
        [~, entry] = min(abs(wave - waveVal(ii)) );
        idx(entry) = 1;
    end
    % Check how we whether the same idx matched two waveVal entries
    nFound = sum(idx);
    if nFound ~= length(waveVal)
        warning('Problems matching wavelengths. Could be out of range.')
    end
end

end