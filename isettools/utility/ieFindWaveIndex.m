function idx = ieFindWaveIndex(wave, waveVal, perfect)
% Returns a (0/1) vector of indices s.t. wave(idx) matches a waveVal entry.
%
% Syntax:
%   idx  = ieFindWaveIndex(wave, waveVal, [perfect])
%
% Description:
%    If we want to address only some wavebands in a radiance data set (e.g., 
%    photons) we find the relevant indices in wave by this call. That is
%      wave(idx)
%    is a vector with the same (or closest match, if perfect == 0)
%    wavelengths.
%
%    Note that wave(idx) is always a column vector, while waveVal could be
%    a row or a column vector.
%
%    If perfect == 1, this routine uses the Matlab function ismember(). 
%
%    If perfect == 0, we accept a closest match, say we want the closest
%    value. In this case, the same wave value may match two waveVal
%    entries, and there will be different vector lengths returned. We
%    announce this mis-match condition.
%
%    Examples contained within the code.
%
% Inputs:
%    wave    - The vector of possible wavelengths.
%    waveVal - The vector of wavelengths whose indices are desired. 
%    perfect - (Optional) Boolean indicating whether or not to look for
%              perfect value match(es) for waveVal in wave. Default is 1. 
%
% Outputs:
%    idx     - A vector of booleans showing which of the instances of wave
%              match the corresponding waveVal entries.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    ieFindWaveIndex

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/30/17  jnm  Formatting & fix example
%    01/16/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    wave = sceneGet(scene, 'wave');
    waveVal = [500, 600]
    idx = ieFindWaveIndex(wave, waveVal);
    wave(idx)
    idx = ieFindWaveIndex(wave, waveVal, false);
    wave(idx)
%}

if notDefined('wave'), error('Must define list of all wavelengths'); end
if notDefined('waveVal'), error('Must define wavelength values'); end
if notDefined('perfect'), perfect = 1; end

if perfect
    % Find only perfect matches
    idx = logical(ismember(wave, waveVal));
else
    % Assume not a member
    idx = false(1, length(wave));   
    
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