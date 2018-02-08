function [idx1, idx2] = ieWave2Index(waveList, wave)
% Convert a wavelength to an index into the wave list.
%
% Syntax:
%   [idx1, idx2] = ieWave2Index(waveList, wave)
%
% Description:
%    If only one return argument is requested, then the index closest to
%    the specified wavelength. If two indices are requested, these are the
%    indices whose wavelength values bound the input wave value. These are
%    always ordered (idx1 < idx2).
%
%    Examples in code.
%
% Inputs:
%    waveList - Wavelength list
%    wave     - Specified wavelength
%
% Outputs:
%    idx1     - If only argument, closest indext to 'wave'. Else, lower
%               bound index around input wave value
%    idx2     - (Optional) If requested, upper bound index around the input
%               wave value
%
% Optional key/value pairs:
%    None.
%
% See also:
%    ieFieldHeight2Index, ieFindWaveIndex
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/21/17  jnm  Formatting
%    01/17/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
	waveList = sceneGet(scene, 'wave');
	idx = ieWave2Index(waveList, 503)
	[idx1, idx2] = ieWave2Index(waveList, 487)
%}

[~, idx1] = min(abs(waveList - wave));

% Determine two indices that bound the wavelength value.
if nargout == 2
    if waveList(idx1) > wave
        % Send back the index below. Order everything properly
        idx2 = max(1, idx1 - 1);
        tmp = idx1;
        idx1 = idx2;
        idx2 = tmp;
    elseif waveList(idx1) < wave
        % Send back the index above. No need to order
        idx2 = min(length(waveList), idx1 + 1);
    else
        idx2 = idx1;
    end
end

end