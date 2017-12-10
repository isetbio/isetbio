function oi = oiSetWvfOptics(oi,wvfP,varargin)
% Put wvf optics into an oi.
%
% Syntax:
%    oi = oiSetWvfOptics(oi,wvf)
%
% Description
%    Take an isetbio wvf struct and put the psf it describes into an oi.
%
%    There is nothing terribly deep here, but this routine takes care of
%    all the fussing.
%
%    This lives in the +ptb package, prepend ptb. when callling. Maybe it
%    should live somewhere else, but it relies on the PTB PsfToOtf
%    function.
%
% Inputs:
%    oi -     Optical image object to update.
%    wvfP -   Wavefront optics structure to insert
%
% Outputs:
%    oi - Updated optical image object.
%
% Optional key/value pairs:
%    'uMPerDegree' -    Scalar, conversion factor between degrees of visual
%                       angle and um on the retina (default 300). You want
%                       this value to match the conversions specified in
%                       the passed oi struct.
%
% See also:
%

% History:
%   12/05/17    dhb   Wrote this.

%% Parse input
p = inputParser;
p.addRequired('oi',@isstruct);
p.addRequired('wvf',@isstruct);
p.addParameter('uMPerDegree', 300, @isscalar);
p.parse(oi,varargin{:});

%% Conversion parameters
uMPerMm = 1000;
uMPerDeg = p.Results.uMPerDegree;

%% Make sure psf is computed on wvf struct and get some parameters we need
wvfP = wvfComputePSF(wvfP);
wvfWls = wvfGet(wvfP, 'calc wave');
nWls = wvfGet(wvfP,'calc nwave');

%% Pull out optics structure
optics = oiGet(oi,'optics');

%% Set optics wavelength
optics = opticsSet(optics,'otf wave',wvfWls);

%% Get otf frequency support and check that it is square
%
% Almost surely true, and haven't thought through all the implications if
% it is not.
sfValuesCyclesMm = opticsGet(optics,'otf support','mm');
if (length(sfValuesCyclesMm{1}) ~= length(sfValuesCyclesMm{2}))
    error('Our code assumes that sf support for otf is square, but it isn''t here.')
end

%% Compute the multispectral otf
for ww = 1:nWls
    
    %
    
    %% Get the gridded spatial frequency support of the otf in cycles/deg.
    %
    % We'll also keep it around in cycles/mm.
    %
    % And convert to support in cycles per degree using 300 um per degree,
    % which is the number that appears to be baked into the optics object.
    
    [xSfGridCyclesMm,ySfGridCyclesMm] = meshgrid(sfValuesCyclesMm{1},sfValuesCyclesMm{2});
    xSfGridCyclesDeg = uMPerDeg*xSfGridCyclesMm/uMPerMm;
    ySfGridCyclesDeg = uMPerDeg*ySfGridCyclesMm/uMPerMm;
    
    %% Get isetbio format OTF at one wavelength
    otf = opticsGet(optics,'otf data',wls(1));
    
    %% Get the desired psf spatial support from the spatial frequency support
    centerPosition = floor(length(sfValuesCyclesMm{1})/2)+1;
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(xSfGridCyclesDeg,ySfGridCyclesDeg);
    position1DMinutes = xGridMinutes(centerPosition,:);
    
    
    
    
    
    % We have to interpolate the psf onto the desired grid in minutes.
    %
    % Set up number of samples for siData PSF and spacing in microns between samples.
    nPSFSamples = p.Results.nPSFSamples;
    umPerSample = p.Results.umPerSample;
    
    iSamp = (1:nPSFSamples) * umPerSample;
    iSamp = iSamp - mean(iSamp);
    iSamp = iSamp(:);
    psf = zeros(nPSFSamples, nPSFSamples, nWls);
    
    %% Do the interplation
    if p.Results.showBar, wBar = waitbar(0, 'Creating PSF'); end
    for ii=1:nWls
        if (p.Results.showBar), waitbar(ii / nWls, wBar); end
        thisPSF = wvfGet(wvfP, 'psf', wls(ii));
        samp = wvfGet(wvfP, 'samples space', 'min', wls(ii));
        samp = samp(:);
        psf(:,:,ii) = interp2(samp, samp', thisPSF, iSamp, iSamp');
    end
    if (p.Results.showBar), close(wBar); end
    
    
    %% Make sure psf has unit volume
    %
    % Not all the routines above guarantee this.
    thePsf = thePsf/sum(thePsf(:));
    
end

%% Stick psf into the optics structure
%
% The ifftshift puts things into the isetbio format.  These are wavelength
% indendent optical estimates.  Not realistic.  We're doing this to compare
% with calculations in the literature that also didn't take chromatic
% aberration into account.
[~,~,theOtfCentered] = PsfToOtf(xGridMinutes,yGridMinutes,thePsf);
theOtfIsetbio = ifftshift(theOtfCentered);
insertOtf = zeros(size(opticsGet(optics,'otf data')));
for ii = 1:length(wls)
    insertOtf(:,:,ii) = theOtfIsetbio;
end
optics = opticsSet(optics,'otf data',insertOtf);

%% Stick optics into oi
oi = oiSet(oi,'optics',optics);
