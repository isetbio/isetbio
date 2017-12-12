function idx = wvfWave2idx(wvf, wList)
% Convert the wavelength list (wList) to indices relative to the list
%
% Syntax:
%   idx = wvfWave2idx(wvf, wList)
%
% Description:
%    Convert the wavelength list (wList) to indices relative to the list
%
%    For example, if wvfGet(wvf, 'wave') is [400 500 600], 
%    and wList is [500, 600] then idx is [2, 3].
%
% Inputs:
%    wvf   - The wavefront object
%    wList - The wavelength list
%
% Outputs:
%    idx   - The indices of wvf that are in the wavelengths
%
% Notes:
%    * [Note: XXX - Currently only returns exact matches within rouding to
%      1 nm. And throws an error if there are none.]
%    * [Note: XXX - Is this really only for 'calc wavelengths'?  Should we
%      have a flag for 'measured wavelength' and 'sce wavelength'?]
%

% History:
%    xx/xx/xx  xxx  wavefront toolbox team
%    11/09/17  jnm  Formatting

% Examples:
%{
  wvf = wvfCreate;
  wvf = wvfSet(wvf, 'wave', 400:10:700);
  wList = 500:100:700;
  idx = wvfWave2idx(wvf, wList)
%}

% Get the wavelengths in the structure
wave = wvfGet(wvf, 'calc wavelengths');

% Check to within 1 nm
idx = find(ismember(round(wave), round(wList)));

% Error if no match
if isempty(idx), error('wvfWave2idx: No matching wavelength in list'); end
return
