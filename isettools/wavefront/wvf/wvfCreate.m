function wvf = wvfCreate(varargin)
% wvf = wvfCreate(varargin)
%
% Create the wavefront parameters structure
%
% varargin:  Structured as param, val pairs
%   
% See also: wvfSet, wvfGet, sceCreate, sceGet
%
% Examples:
%    wvf = wvfCreate('wave',[400:10:700]);
%
% History
%   7/20/12 dhb      Get rid of weighting spectrum, replace with cone psf info structure
%
% (c) Wavefront Toolbox Team 2011, 2012

%% Book-keeping
wvf = [];
wvf = wvfSet(wvf,'name','default');
wvf = wvfSet(wvf,'type','wvf');

%% Zernike coefficients and related
wvf = wvfSet(wvf,'zcoeffs',0);
wvf = wvfSet(wvf,'measured pupil',8);
wvf = wvfSet(wvf,'measured wl',550);
wvf = wvfSet(wvf,'measured optical axis',0);
wvf = wvfSet(wvf,'measured observer accommodation',0);
wvf = wvfSet(wvf,'measured observer focus correction',0);

%% Spatial sampling parameters
wvf = wvfSet(wvf,'sample interval domain','psf');
wvf = wvfSet(wvf,'spatial samples',201);
wvf = wvfSet(wvf,'ref pupil plane size',16.212);

%% Calculation parameters
wvf = wvfSet(wvf,'calc pupil size',3);
wvf = wvfSet(wvf,'calc wavelengths',550);
wvf = wvfSet(wvf,'calc optical axis',0);
wvf = wvfSet(wvf,'calc observer accommodation',0);
wvf = wvfSet(wvf,'calc observer focus correction',0);

% Cone sensitivities and weighting spectrum for combining the PSFs across
% wavelengths.  We keep these as a structure at something resembling a
% wide range of wavelength samples, along with the wavelength info.  
% When we use them, we spline down to the wavelength sampling of the psfs.
% These follow PTB spectral conventions.
load('T_cones_ss2');
conePsfInfo.S = S_cones_ss2;
conePsfInfo.T = T_cones_ss2;
conePsfInfo.spdWeighting = ones(conePsfInfo.S(3),1);
wvf = wvfSet(wvf,'calc cone psf info',conePsfInfo);

% Stiles Crawford Effect parameters. 
wvf = wvfSet(wvf,'sce params',sceCreate([],'none'));

% Handle any additional arguments via wvfSet
if ~isempty(varargin)
    if isodd(length(varargin))
        error('Arguments must be (pair, val) pairs');
    end
    for ii=1:2:(length(varargin)-1)
        wvf = wvfSet(wvf,varargin{ii},varargin{ii+1});
    end
end

return
