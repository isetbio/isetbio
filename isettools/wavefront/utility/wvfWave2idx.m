function idx = wvfwvfWave2idx(wvf,wList)
% idx = wvfwvfWave2idx(wvf,wList)
%
% Convert the wavelength list (wList) to indices relative to
%
% For example, if wvfGet(wvf,'wave') is [400 500 600],
% and wList is [500,600] then idx is [2,3].
%
% Note: Currently only returns exact matches within rouding to 1 nm.
% And throws an error if there are none.
%
% Note: Is this really only for 'calc wavelengths'?  Should we have a flag for
% 'measured wavelength' and 'sce wavelength'?
%
% Example
%  wvf = wvfCreate;
%  wvf = wvfSet(wvf,'wave',400:10:700);
%  wList = 500:100:700;
%  idx = wvfwvfWave2idx(wvf,wList)

% Get the wavelengths in the structure
wave = wvfGet(wvf,'calc wavelengths');

% Check to within 1 nm
idx = find(ismember(round(wave),round(wList)));

% Error if no match
if isempty(idx), error('wvfWave2idx: No matching wavelength in list'); end
return
