% t_opticsGetAndSetPsf
%
% Description:
%   First, show how to get and set psf's in units of minutes of arc, or from otfs
%   specified over frequencies in cycles/deg. This ability is useful when
%   we want to build isetbio optics using various estimates of optical quality
%   in the literature.
%
%   Then show how to do this using the wavefront optics.
%
%   This tutorial also serves some validation purposes, as it checks that
%   various ways of doing the same thing give the same answer.
%
% See also:
%

% History:
%   2/5/17  dhb  Wrote starting with an extant tutorial.

%% Initialize
clear; ieInit;

%% Create the oi structure and pull out the optics
%
% Snag the wavelengths while we're at it.
theWl = 550;
oi = oiCreate('wvf human');
optics = oiGet(oi,'optics');
wls = opticsGet(optics,'wave');

%% Check that support is square
% 
% Almost surely true, and haven't thought through all the implications if
% it is not.
sfValuesCyclesMm = opticsGet(optics,'otf support','mm');
if (length(sfValuesCyclesMm{1}) ~= length(sfValuesCyclesMm{2}))
    error('Our code assumes that sf support for otf is square, but it isn''t here.')
end

%% Change even/odd support to odd/even support
%
% If dimension of otf/psf support is odd, make it even.  If it is even, 
% make it odd.  This is a bit of a hack but lets us test that things work
% for both odd and even dimension.
%
% The case where it starts off even is not yet handled.
CHANGE_SUPPORTDIM = false;
if (CHANGE_SUPPORTDIM)
    % Case where it starts odd
    if (rem(length(sfValuesCyclesMm{1}),2) ~= 0)
        
        % Lop off highest postive frequency in both x and y and put back.
        sfValuesCyclesMm{1} = sfValuesCyclesMm{1}(1:end-1);
        optics = opticsSet(optics,'otffx',sfValuesCyclesMm{1});
        sfValuesCyclesMm{2} = sfValuesCyclesMm{2}(1:end-1);
        optics = opticsSet(optics,'otffy',sfValuesCyclesMm{2});
        
        % If it starts out odd, lop off highest postive frequency in the otf.
        % Isetbio stores the otf in Matlab's first entry is zero freuqency
        % format, but the lopping is easier to think about in the zero
        % frequency at the center format. Use fftshift and ifftshift to go back
        % and forth, and do the lopping in between.  We need to do all the
        % wavelengths, and put the whole resized cube back in, in one fell
        % swoop.
        otfOddSupport = opticsGet(optics,'otf data');
        for ii = 1:length(wls)
            otfCentered = fftshift(otfOddSupport(:,:,ii));
            otfCentered = otfCentered(1:end-1,1:end-1);
            otfEvenSupport(:,:,ii) = ifftshift(otfCentered);
        end
        optics = opticsSet(optics,'otf data',otfEvenSupport);
        clear otfOddSupport otfEvenSupport otfCentered
    % Case where it starts even
    else
        error('Not yet implemented.');
    end
end

%% Get the gridded spatial frequency support of the otf in cycles/deg.
%
% We'll also keep it around in cycles/mm.
%
% And convert to support in cycles per degree using 300 um per degree,
% which is the number that appears to be baked into the optics object.
uMPerMm = 1000;
uMPerDegree = 300;
[xSfGridCyclesMm,ySfGridCyclesMm] = meshgrid(sfValuesCyclesMm{1},sfValuesCyclesMm{2});
xSfGridCyclesDegree = uMPerDegree*xSfGridCyclesMm/uMPerMm;
ySfGridCyclesDegree = uMPerDegree*ySfGridCyclesMm/uMPerMm;

%% Get isetbio format OTF at a specified wavelength
otf = opticsGet(optics,'otf data',theWl);

%% Derive the psf from the otf
%
% We have to convert to the zero sf at center representation to use
% OtfToPsf, using fftshift.
[xGridMinutes,yGridMinutes,psf] = OtfToPsf(xSfGridCyclesDegree,ySfGridCyclesDegree,fftshift(otf));
centerPosition = floor(length(sfValuesCyclesMm{1})/2)+1;
position1DMinutes = xGridMinutes(centerPosition,:);
wvfHuman1DPsf = psf(centerPosition,:);

%% Show otf and psf in 2D
figure; clf;
subplot(1,2,1); hold on;
mesh(xSfGridCyclesDegree,ySfGridCyclesDegree,fftshift(abs(otf)));
xlim([-100 100]); ylim([-100 100]);
axis('square');
xlabel('X SF (cycles/deg)');
ylabel('Y SF (cycles/deg)');
title('OTF');
subplot(1,2,2); hold on
mesh(xGridMinutes,yGridMinutes,psf);
xlim([-10 10]); ylim([-10 10]);
axis('square');
xlabel('x (minutes)');
ylabel('y (minutes)');
title('PSF');

%% Get Davila-Geisler PSF and plot against what isetbio wvf looks like
DavilaGeislerLsf = DavilaGeislerLSFMinutes(position1DMinutes);
DavilaGeislerPsf = LsfToPsf(DavilaGeislerLsf);
psfFig = figure; clf; hold on
plot(position1DMinutes,wvfHuman1DPsf/max(wvfHuman1DPsf),'r','LineWidth',4);
plot(position1DMinutes,DavilaGeislerPsf(centerPosition,:)/max(DavilaGeislerPsf(centerPosition,:)),'g-','LineWidth',4);
xlim([-4 4]);
xlabel('Position (minutes');
ylabel('Normalized PSF Slice');
title('PSF')
legend({sprintf('Wvf Human @%d nm',theWl),'Davila-Geisler'});

%% Stick DavilaGeisler into the optics structure
%
% The ifftshift puts things into the isetbio format.
[~,~,DavilaGeislerOtfCentered] = PsfToOtf(xGridMinutes,yGridMinutes,DavilaGeislerPsf);
DavilaGeislerOtfIsetbio = ifftshift(DavilaGeislerOtfCentered);
insertOtf = zeros(size(opticsGet(optics,'otf data')));
for ii = 1:length(wls)
    insertOtf(:,:,ii) = DavilaGeislerOtfIsetbio;
end
optics = opticsSet(optics,'otf data',insertOtf);

%% Make sure everything is hunky-dory by making the plot using isetbio's fcn
oi = oiSet(oi,'optics',optics);
udata = oiPlot(oi,'psf',[],theWl);
figure(psfFig);
plot(60*udata.x(centerPosition,:)/uMPerDegree,udata.psf(centerPosition,:)/max(udata.psf(centerPosition,:)),'k-','LineWidth',2);

%% Let's do this through the routine and show that we get the same answer
oi2 = oiCreate('wvf human');
oi2 = ptb.oiSetPtbOptics(oi2);
udata2 = oiPlot(oi2,'psf',[],theWl);
figure(psfFig);
plot(60*udata2.x(centerPosition,:)/uMPerDegree,udata2.psf(centerPosition,:)/max(udata2.psf(centerPosition,:)),'r:','LineWidth',2);

%% Add Westheimer and Williams estimates for fun
oi2 = ptb.oiSetPtbOptics(oi2,'opticsModel','Westheimer');
udata2 = oiPlot(oi2,'psf',[],theWl);
figure(psfFig);
plot(60*udata2.x(centerPosition,:)/uMPerDegree,udata2.psf(centerPosition,:)/max(udata2.psf(centerPosition,:)),'c','LineWidth',2);
oi2 = ptb.oiSetPtbOptics(oi2,'opticsModel','Williams');
udata2 = oiPlot(oi2,'psf',[],theWl);
figure(psfFig);
plot(60*udata2.x(centerPosition,:)/uMPerDegree,udata2.psf(centerPosition,:)/max(udata2.psf(centerPosition,:)),'y','LineWidth',2);

%% And finally plot the Davila-Geisler lsf if we take it directly as the psf.
oi2 = ptb.oiSetPtbOptics(oi2,'opticsModel','DavilaGeislerLsfAsPsf');
udata2 = oiPlot(oi2,'psf',[],theWl);
figure(psfFig);
plot(60*udata2.x(centerPosition,:)/uMPerDegree,udata2.psf(centerPosition,:)/max(udata2.psf(centerPosition,:)),'k:','LineWidth',2);
legend({sprintf('Wvf Human @%d nm',theWl),'Davila-Geisler','D/G Again','D/G Yet Again','Westheimer','Williams','D/G Lsf as Psf'});

%% Now start again, and work with wavefront optics PSFs
%
% Use Thibos measurements, but set some defocus to make
% the PSF more interesting.
%
% To get all the ways to come out consistently, we need to be careful
% to use consistent parameters across all of them, so we take some care
% to define and set spatial sampling for the PSF here, as well as
% umPerDegree.
pupilMM = 6; zCoeffs = wvfLoadThibosVirtualEyes(pupilMM);
psfSpatialSamples = 200;
psfUmPerSample = 0.25;
umPerDegree = 300;
psfMinPerSample = 60*psfUmPerSample/umPerDegree;
wvfP = wvfCreate('calc wavelengths', theWl, ...
        'zcoeffs', zCoeffs, 'measured pupil', pupilMM, ...
        'spatialSamples',psfSpatialSamples, ...
        'umPerDegree',umPerDegree, ...
        'name', sprintf('human-%d', pupilMM));
wvfP = wvfSet(wvfP,'ref psf sample interval',psfMinPerSample);
wvfP = wvfSet(wvfP,'zcoeffs',1,'defocus');
wvfP = wvfComputePSF(wvfP);

% Get and plot the psf obtained directly from the wvf structure.
% This is what we ought to get back from an oi/optics structure,
% if we put it in correctly.
psfFig3 = figure; hold on
psf3FromWvf = wvfGet(wvfP,'1d psf',theWl);
psf3FromWvfSpatialSamples1D = wvfGet(wvfP,'psf angular samples','min',theWl);
plot(psf3FromWvfSpatialSamples1D, ...
    psf3FromWvf/max(psf3FromWvf),...
    'c','LineWidth',6);

% Convert to oi using wvf2oi.  When we get the psf data
% using oiPlot, the support is in microns.  Convert to 
% minutes and plot.
oi3 = wvf2oi(wvfP);
udata3 = oiPlot(oi3,'psf',[],theWl);
supportRowSize = size(udata3.x,1);
centerPosition = floor(supportRowSize/2)+1;
plot(60*udata3.x(centerPosition,:)/uMPerDegree, ...
    udata3.psf(centerPosition,:)/max(udata3.psf(centerPosition,:)),...
    'b','LineWidth',4);

% Get isetbio format OTF back out of the oi struct, at the specified wavelength
optics3 = oiGet(oi3,'optics');
otf3 = opticsGet(optics3,'otf data',theWl);

% Derive the psf back from the otf using the PTB routine OtfToPsf.  As
% often happens, we have independently done the same things in isetbio and
% PTB, and here we want to make sure that we get the same answer.
%
% Before calling the PTB routine OtfToPsf, we have to convert to the zero
% sf at center representation, using fftshift.
sfValuesCyclesMm3 = opticsGet(optics3,'otf support','mm');
[xSfGridCyclesMm3,ySfGridCyclesMm3] = meshgrid(sfValuesCyclesMm3{1},sfValuesCyclesMm3{2});
xSfGridCyclesDegree3 = uMPerDegree*xSfGridCyclesMm3/uMPerMm;
ySfGridCyclesDegree3 = uMPerDegree*ySfGridCyclesMm3/uMPerMm;
[xGridMinutes3,yGridMinutes3,psf3] = OtfToPsf(xSfGridCyclesDegree3,ySfGridCyclesDegree3,fftshift(otf3));
centerPosition3 = floor(length(sfValuesCyclesMm3{1})/2)+1;
position1DMinutes3 = xGridMinutes3(centerPosition3,:);
wvfHuman1DPsf3 = psf3(centerPosition,:);
figure(psfFig3);
plot(position1DMinutes3,wvfHuman1DPsf3/max(wvfHuman1DPsf3),'r','LineWidth',2);

%% Not yet working
% % Do the conversion using si format
% [siPSFData, wvfP] = wvf2PSF(wvfP,'nPSFSamples',psfSpatialSamples, ...
%     'umPerSample',psfUmPerSample, ...
%     'showBar',false);
% 
% % Convert to optics and then oi using siSynthetic
% oi4 = oiCreate('human');
% optics4 = oiGet(oi4,'optics');
% optics4 = opticsSet(optics4,'otf',zeros(psfSpatialSamples,psfSpatialSamples));
% oi4 = oiSet(oi4,'optics',optics4);
% optics4 = siSynthetic('custom', oi4, siPSFData);
% oi4 = oiSet(oi4, 'optics', optics4);
% 
% udata4 = oiPlot(oi4,'psf',[],theWl);
% supportRowSize4 = size(udata4.x,1);
% centerPosition4 = floor(supportRowSize4/2)+1;
% figure(psfFig3);
% plot(60*udata4.x(centerPosition4,:)/uMPerDegree, ...
%     udata4.psf(centerPosition4,:)/max(udata4.psf(centerPosition4,:)),...
%     'k','LineWidth',1);
% xlim([-15 15]);
% 
% % flength = 0.017;  % Human focal length is 17 mm
% % oi = oiSet(oi, 'optics fnumber', flength/pupilMM);
% % oi = oiSet(oi, 'optics flength', flength);



